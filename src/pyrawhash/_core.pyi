from __future__ import annotations
import typing

__all__ = ['Alignment', 'Index', 'Request', 'Response', 'ResponseGenerator']

class Alignment:
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, header: str, fwd: bool, start: typing.SupportsInt, pres_frac: typing.SupportsFloat, qry_len: typing.SupportsInt) -> None:
        ...
    @property
    def ctg(self) -> str:
        ...
    @property
    def pres_frac(self) -> float:
        ...
    @property
    def r_en(self) -> int:
        ...
    @property
    def r_st(self) -> int:
        ...
    @property
    def strand(self) -> int:
        ...
class Index:
    def __init__(self, argv: list) -> None:
        """
        Load a spumoni index.
        :param argv: commandline arguments
        """
        ...
    def validate(self): -> None:
    ...
def query(self, sequence: str) -> Alignment:
    """
    Query a sequence in the index
    :param sequence: query sequence
    :return: an alignment of the query
    """
    ...
def query_stream(self, requests: typing.Iterator) -> ResponseGenerator:
    """
    Query a stream of requests and return a stream of responses (Readfish compatible)
    :param requests: an iterator of query request objects
    :return: A generator of query response objects
    """
    ...