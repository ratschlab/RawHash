from __future__ import annotations
from typing import Iterable, TYPE_CHECKING

from collections import namedtuple

import numpy as np
import numpy.typing as npt
import minknow_api
from minknow_api.protocol_pb2 import ProtocolRunInfo
from minknow_api.device_pb2 import GetSampleRateResponse

from readfish.plugins.utils import Result

CALIBRATION = namedtuple("calibration", "scaling offset")

class DefaultDAQValues:
    """Provides default calibration values

    Mimics the read_until_api calibration dict value from
    https://github.com/nanoporetech/read_until_api/blob/2319bbe/read_until/base.py#L34
    all keys return scaling=1.0 and offset=0.0
    """

    calibration = CALIBRATION(1.0, 0.0)

    def __getitem__(self, _):
        return self.calibration


_DefaultDAQValues = DefaultDAQValues()

class Caller:
    def __init__(self, debug_log: str | None = None, **kwargs) -> None:
        print("No-op basecaller initialized")

    def validate(self) -> None:
        pass

    def basecall(self, reads: list[tuple[int, minknow_api.data_pb2.GetLiveReadsResponse.ReadData]],
                 signal_dtype: np.dtype, daq_values: dict[int, CALIBRATION]) -> Iterable[Result]:
        daq_values = _DefaultDAQValues if daq_values is None else daq_values
        for channel, read in reads:
            raw_data=np.frombuffer(read.raw_data, signal_dtype)
            daq_offset=daq_values[channel].offset,
            daq_scaling=daq_values[channel].scaling,
            if signal_dtype.kind != 'f':
                raw_data = (raw_data + daq_offset) * daq_scaling
            yield Result(channel=channel, read_id=read.id, seq='', basecall_data=raw_data)

    def describe(self) -> str:
        return "This caller does nothing except appending the raw signal to Result.basecall_data"

    def disconnect(self) -> None:
        pass