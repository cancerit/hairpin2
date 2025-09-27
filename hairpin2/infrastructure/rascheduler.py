from collections.abc import Iterable
from typing import Any, Self, cast

from hairpin2.infrastructure.process import ReadAwareProcess, ReadAwareProcessProtocol
from hairpin2.infrastructure.process_engines import EngineResult_T, FlagResult_T, ProcessTypeEnum
from hairpin2.infrastructure.process_params import RunParams
from hairpin2.infrastructure.structures import FlagResult


class RAExec:
    _mandate_excludes: bool
    _taggers: tuple[ReadAwareProcessProtocol[None], ...]
    _flaggers: tuple[ReadAwareProcessProtocol[FlagResult], ...]

    def __init__(
        self,
        taggers_in_order: Iterable[ReadAwareProcessProtocol[None]],
        flaggers_in_order: Iterable[ReadAwareProcessProtocol[FlagResult_T]],
        mandate_excludes: bool,
    ):
        self._mandate_excludes = mandate_excludes
        _, tags_added = self.validate_tagger_order(
            taggers_in_order, mandate_excludes, raise_on_fail=True
        )
        _ = self.validate_flagger_exec(
            flaggers_in_order, tags_added, mandate_excludes, raise_on_fail=True
        )
        self._taggers = tuple(taggers_in_order)
        self._flaggers = tuple(flaggers_in_order)

    @property
    def mandate_excludes(self):
        return self._mandate_excludes

    @property
    def taggers(self):
        return self._taggers

    @property
    def flaggers(self):
        return self._flaggers

    @staticmethod
    def _partition_processes(
        processes: Iterable[ReadAwareProcessProtocol[EngineResult_T]], raise_on_fail: bool = False
    ) -> tuple[
        bool,
        tuple[
            tuple[ReadAwareProcessProtocol[None], ...],
            tuple[ReadAwareProcessProtocol[FlagResult], ...],
        ],
    ]:
        flaggerl: list[ReadAwareProcessProtocol[FlagResult]] = []
        taggerl: list[ReadAwareProcessProtocol[None]] = []
        taggers_all_registered = False
        valid = True
        for proc in processes:
            match proc.ProcessType:
                case ProcessTypeEnum.TAGGER:
                    if taggers_all_registered:
                        if raise_on_fail:
                            raise RuntimeError("tagger appeared in execution order after flagger")
                        valid = False
                    taggerl.append(cast(ReadAwareProcessProtocol[None], proc))
                case ProcessTypeEnum.FLAGGER:
                    if not taggers_all_registered:
                        taggers_all_registered = True
                    flaggerl.append(cast(ReadAwareProcessProtocol[FlagResult], proc))
                case _:
                    raise TypeError  # pyright: ignore[reportUnreachable]
        return valid, (tuple(taggerl), tuple(flaggerl))

    @staticmethod
    def check_namespacing_clashfree(
        processes: Iterable[ReadAwareProcessProtocol[Any]]
        | Iterable[type[ReadAwareProcessProtocol[Any]]],
    ):
        nsl = [proc.ProcessNamespace for proc in processes]
        return len(nsl) == len(set(nsl))  # clashes if false

    @staticmethod
    def validate_exec(
        processes: Iterable[ReadAwareProcessProtocol[EngineResult_T]],
        mandate_excludes: bool,
        raise_on_fail: bool,
    ) -> tuple[
        bool,
        tuple[
            tuple[ReadAwareProcessProtocol[None], ...],
            tuple[ReadAwareProcessProtocol[FlagResult], ...],
        ],
    ]:
        valid, part_procs = RAExec._partition_processes(processes, raise_on_fail)
        if valid:
            valid, tags_added = RAExec.validate_tagger_order(
                part_procs[0], mandate_excludes, raise_on_fail
            )
            if valid:
                valid = RAExec.validate_flagger_exec(
                    part_procs[1], tags_added, mandate_excludes, raise_on_fail
                )
        return valid, part_procs

    # make wrapped read register itself as marked in it's parent readview if there is one
    @staticmethod
    def validate_tagger_order(
        processes: Iterable[ReadAwareProcessProtocol[None]],
        mandate_excludes: bool,
        raise_on_fail: bool = False,
    ) -> tuple[bool, set[str]]:
        """
        validate that the order of tagger processes as provided in the processses iterable can be run without error.
        In pratice this means checking that at a minimum, all marks/tags listed in a given process as required (i.e.
        in require_marks) have been checked and applied by processes occuring prior in the run order.

        Args:
            processes (Iterable): Iterable of processes, the order of which defines the run order
            mandate_excludes (bool): if true, then for each process, all marks listed as exclude_marks must have been applied by tagger processes coming prior in the run order, in addition to require_marks
            raise_on_fail (bool): whether to throw an error on recieving an illegal order, or just return a bool
        """
        tags_set: set[str] = set()
        valid = True
        for proc in processes:
            assert proc.AddsMarks is not None
            reqs_fulfilled = proc.require_marks.issubset(tags_set)
            if mandate_excludes:
                reqs_fulfilled = reqs_fulfilled and proc.exclude_marks.issubset(tags_set)
            if not reqs_fulfilled:
                if raise_on_fail:
                    raise RuntimeError("Invalid tagger execution order")
                valid = False
                break
            if tags_set & proc.AddsMarks:  # if overlap
                if raise_on_fail:
                    raise RuntimeError("taggers appear to add same tag. This is not supported")
                valid = False
                break
            tags_set = tags_set | proc.AddsMarks
        return (valid, tags_set) if valid else (valid, set())

    @staticmethod
    def validate_flagger_exec(
        processes: Iterable[ReadAwareProcessProtocol[FlagResult]],
        tags_set: set[str],
        mandate_excludes: bool,
        raise_on_fail: bool = False,
    ):
        valid = True
        for proc in processes:
            reqs_fulfilled = proc.require_marks.issubset(tags_set)
            if mandate_excludes:
                reqs_fulfilled = reqs_fulfilled and proc.exclude_marks.issubset(tags_set)
            if not reqs_fulfilled:
                if raise_on_fail:
                    raise RuntimeError("Invalid tagger execution order")
                valid = False
                break
        return valid

    @classmethod
    def from_config(
        cls,
        configd: dict[str, Any],
        proc_pool: Iterable[type[ReadAwareProcess]],
    ) -> Self:
        """
        initialise read-aware executor from config dict.

        currently fragile

        Args:
            configd (dict[str, Any]): Dictionary describing how processes should execute.
            proc_pool (Iterable[type[ReadAwareProcess]]): Iterable of process definitions which may be instatiated for a given run.
        """
        if not all(isinstance(proc, ReadAwareProcessProtocol) for proc in proc_pool):
            raise ValueError(
                "One or more processes provided to proc_pool do not have the appropriate methods"
            )
        cast_proc_pool = cast(Iterable[type[ReadAwareProcessProtocol[Any]]], proc_pool)

        const_config_keys = ["params", "exec"]
        if not set(const_config_keys) == configd.keys():
            raise ValueError(
                f"config does not have correct top-level keys, expected {const_config_keys}, recieved {list(configd.keys())}"
            )

        excludes_opt = cast(bool, configd["exec"]["opts"]["mandate-excludes"])
        paramd = cast(dict[str, Any], configd["params"])
        execd = cast(dict[str, Any], configd["exec"])
        if not isinstance(excludes_opt, bool):
            raise ValueError
        if not isinstance(execd, dict):
            raise ValueError
        if not isinstance(paramd, dict):
            raise ValueError

        if not cls.check_namespacing_clashfree(cast_proc_pool):
            raise ValueError("Processes with clashing namespaces")

        proc_nsd = {proc_type.ProcessNamespace: proc_type for proc_type in cast_proc_pool}

        if not set(paramd.keys()).issubset(proc_nsd):
            raise ValueError("config requesting unknown processes")

        # TODO: if using pydantic doesn't engender end user having to necessarily define paramater dataclasses,
        # then clearly need pydantic validated config dataclass (may suffice in the meantime anyway since this is hp2 not htsflow)
        # config order defines execution order
        ordered_proc_inst = tuple(
            proc_nsd[proc_name](
                paramd.get(proc_name, {}),
                execd[proc_name]["require-marks"],
                execd[proc_name]["exclude-marks"],
            )
            for proc_name in execd.keys()
            if proc_name != "opts" and execd[proc_name]["enable"]
        )

        _, (taggers, flaggers) = cls.validate_exec(
            ordered_proc_inst, excludes_opt, raise_on_fail=True
        )

        # this duplicates some of the above validation at the moment. Oh well.
        ra_ex = cls(taggers, flaggers, excludes_opt)
        return ra_ex

    def run(self, run_data: RunParams) -> tuple[FlagResult, ...]:
        for proc in self.taggers:
            proc(run_data)
            proc.reset()

        res: list[FlagResult] = []
        for proc in self.flaggers:
            res.append(proc(run_data))
            proc.reset()

        return tuple(res)
