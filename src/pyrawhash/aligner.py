import os
import sys
from pathlib import Path
from typing import Optional, Any, Iterable

import attrs

from . import Index, Request, Alignment

def convert_to_argv(*args, **kwargs):
    """
    Converts Python *args and **kwargs into a sys.argv-style list.

    This function is designed to build a list of command-line arguments
    from Python function arguments.

    - `kwargs` are converted into long-form options (e.g., `threads=8` -> `['--threads', '8']`).
    - Underscores (`_`) in kwarg keys are converted to hyphens (`-`) (e.g., `doc_array=True` -> `['--doc-array']`).
    - Boolean `True` values add only the flag (e.g., `MS=True` -> `['--MS']`).
    - Boolean `False` values cause the flag to be omitted entirely.
    - All other values are converted to strings and appended after their flag.
    - `args` are appended as positional arguments at the end of the list.

    Returns:
        list: A list of strings formatted like sys.argv (minus the script name).
    """
    argv_list = ['pyrawhash']

    # --- Handle Keyword Arguments (Flags & Options) ---
    for key, value in kwargs.items():
        if len(key) == 1:
            flag = '-' + key        # e.g. v=True  -> '-v'
        else:
            # 'small_window' -> '--small-window'
            flag = '--' + key.replace('_', '-')

        if value is True:
            # Handle boolean flags like --MS or --classify
            argv_list.append(flag)
        elif value is False:
            # If a flag is explicitly False, skip it.
            # This allows setting defaults to False.
            pass
        else:
            # Handle options that take a value, like --threads 8
            argv_list.append(flag)
            argv_list.append(str(value))

    # --- Handle Positional Arguments ---
    # These are typically input files or other arguments
    # that don't have flags.
    for arg in args:
        argv_list.append(str(arg))

    print("Argv:", argv_list)
    return argv_list

@attrs.define
class Result:
    """Result holder

    This should be progressively filled with data from the basecaller,
    barcoder, and then the aligner.

    :param channel: The channel that this read is being sequenced on
    :param read_id: The read ID assigned to this read by MinKNOW
    :param seq: The basecalled sequence for this read
    :param decision: The ``Decision`` that has been made, this will by used to determine the ``Action``
    :param barcode: The barcode that has been assigned to this read
    :param basecall_data: Any extra data that the basecaller may want to send to the aligner
    :param alignment_data: Any extra alignment data
    """

    channel: int
    read_id: str
    seq: str
    barcode: Optional[str] = attrs.field(default=None)
    basecall_data: Optional[Any] = attrs.field(default=None)
    alignment_data: Optional[list[Alignment]] = attrs.field(default=None)

class Aligner:

    def __init__(self, debug_log: str | None = None, **kwargs):
        if debug_log:
            if debug_log == 'stdout':
                self.logfile = sys.stdout
            elif debug_log == 'stderr':
                self.logfile = sys.stderr
            else:
                self.logfile = open(debug_log, 'w')
        else:
            self.logfile = None
        if 'threads' in kwargs:
            n_threads = int(kwargs['threads'])
            if n_threads > 0:
                os.environ['PARLAY_NUM_THREADS'] = str(n_threads)
            kwargs.pop('threads')
        else:
            os.environ['PARLAY_NUM_THREADS'] = '1'

        if 'idx' not in kwargs:
            print("Must provide an index using the '--idx' parameter")
            exit(1)
        else:
            self.idx = kwargs['idx']
            kwargs.pop('idx')

        self.kwargs = kwargs
        args_list = convert_to_argv(**kwargs)
        args_list.append(self.idx)
        print("Args:", args_list)
        self.aligner = Index(args_list)

    def validate(self) -> None:
        self.aligner.validate()
        print("PyRawhash: All OK!")

    @property
    def initialised(self) -> bool:
        return True

    def describe(self, regions: list, barcodes: dict) -> str:
        return ("Loaded index from prefix {}. Hopefully, I will return a more meaningful description in the future".
                format(self.idx))

    def map_reads(self, calls: Iterable[Result]) -> Iterable[Result]:
        skipped = []
        metadata = {}
        def _gen(_basecalls):
            """Create request objects for aligner."""
            for result in _basecalls:
                id = result.read_id
                metadata[id] = result
                seq = result.seq
                if not seq:
                    skipped.append(result)
                    continue
                yield Request(channel=result.channel, id=id, signal=result.basecall_data)

        responses = self.aligner.query_stream(_gen(calls))
        for response in responses:
            result = metadata[response.id]
            result.alignment_data = [response.alignment] if response.alignment.ctg != '*' else []
            yield result
        for result in skipped:
            result.alignment_data = []
            yield result

    def disconnect(self):
        if self.logfile and self.logfile not in [sys.stdout, sys.stderr]:
            self.logfile.close()