from collections.abc import Iterable
from typing import Any, Self, cast

from hairpin2.abstractions.process import ReadAwareProcess
from hairpin2.abstractions.process_engines import ProcessTypeEnum
from hairpin2.abstractions.process_params import RunParams
# pyright: reportAny=false
# pyright: reportExplicitAny=false
# pyright: reportUnnecessaryIsInstance=false

class RAExec:
    _mandate_excludes: bool
    _taggers: tuple[ReadAwareProcess, ...]
    _flaggers: tuple[ReadAwareProcess, ...]

    def __init__(
        self,
        register_processes: Iterable[ReadAwareProcess] | None,
        mandate_excludes: bool
    ):
        self._mandate_excludes = mandate_excludes
        if register_processes is not None:
            _, part_procs = RAExec.validate_exec(register_processes, mandate_excludes, raise_on_fail=True)

            # next figure out which is most important, schedule that first
            # idk reduce into dependency graph of processes
            # OK DON'T DO THIS
            # as Luca points out, at the very least as a first pass,
            # just figure out if the order as input is valid and fail if not
            # also requires/excludes can be defined in the config
            self._taggers = part_procs[0]
            self._flaggers = part_procs[1]
        else:
            self._taggers = ()
            self._flaggers = ()

    @property
    def mandate_excludes(
        self
    ):
        return self._mandate_excludes

    @property
    def taggers(
        self
    ):
        return self._taggers

    @property
    def flaggers(
        self
    ):
        return self._flaggers

    @staticmethod
    def _partition_processes(
        processes: Iterable[ReadAwareProcess],
        raise_on_fail: bool = False
    ) -> tuple[bool, tuple[tuple[ReadAwareProcess, ...], tuple[ReadAwareProcess, ...]]]:
        flaggerl: list[ReadAwareProcess] = []
        taggerl: list[ReadAwareProcess] = []
        taggers_all_registered = False
        valid = True
        for proc in processes:
            match proc.ProcessType:
                case ProcessTypeEnum.TAGGER:
                    if taggers_all_registered:
                        if raise_on_fail:
                            raise RuntimeError('tagger appeared in execution order after flagger')
                        valid = False
                    taggerl.append(proc)
                case ProcessTypeEnum.FLAGGER:
                    if not taggers_all_registered:
                        taggers_all_registered = True
                    flaggerl.append(proc)
                case _:
                    raise TypeError
        return valid, (tuple(taggerl), tuple(flaggerl))

    @staticmethod
    def check_namespacing_clashfree(
        processes: Iterable[ReadAwareProcess] | Iterable[type[ReadAwareProcess]]
    ):
        nsl = [proc.ProcessNamespace for proc in processes]
        return len(nsl) == len(set(nsl))  # clashes if false

    @staticmethod
    def validate_exec(
        processes: Iterable[ReadAwareProcess],
        mandate_excludes: bool,
        raise_on_fail: bool
    ) -> tuple[bool, tuple[tuple[ReadAwareProcess, ...], tuple[ReadAwareProcess, ...]]]:
        valid, part_procs = RAExec._partition_processes(processes, raise_on_fail)
        if valid:
            valid, tags_added = RAExec.validate_tagger_order(part_procs[0], mandate_excludes, raise_on_fail)
            if valid:
                valid = RAExec.validate_flagger_exec(part_procs[1], tags_added, mandate_excludes, raise_on_fail)
        return valid, part_procs

    # need to make clearer pattern of boolean read marks
    # and data holding read marks
    # and make wrapped read register itself as marked in it's parent readview if there is one
    @staticmethod
    def validate_tagger_order(
        processes: Iterable[ReadAwareProcess],
        mandate_excludes: bool,
        raise_on_fail: bool = False
    ) -> tuple[bool, set[str]]:
        tags_set: set[str] = set()
        valid = True
        for proc in processes:
            reqs_fulfilled = proc.require_marks.issubset(tags_set)
            if mandate_excludes:
                reqs_fulfilled = reqs_fulfilled and proc.exclude_marks.issubset(tags_set)
            if not reqs_fulfilled:
                if raise_on_fail:
                    raise RuntimeError('Invalid tagger execution order')
                valid = False
                break
            if tags_set & proc.AddsMarks:  # if overlap  # TODO/BUG: AddsMarks might have to be a var on ReadAwareProcess... or monkey patched onto the engine. Or go back to init subclass, but long
                if raise_on_fail:
                    raise RuntimeError('taggers appear to add same tag. This is not supported')
                valid = False
                break
            tags_set = tags_set | proc.AddsMarks
        return (valid, tags_set) if valid else (valid, set())

    @staticmethod
    def validate_flagger_exec(
        processes: Iterable[ReadAwareProcess],
        tags_set: set[str],
        mandate_excludes: bool,
        raise_on_fail: bool = False
    ):
        valid = True
        for proc in processes:
            reqs_fulfilled = proc.require_marks.issubset(tags_set)
            if mandate_excludes:
                reqs_fulfilled = reqs_fulfilled and proc.exclude_marks.issubset(tags_set)
            if not reqs_fulfilled:
                if raise_on_fail:
                    raise RuntimeError('Invalid tagger execution order')
                valid = False
                break
        return valid

    @classmethod
    def from_config(
        cls,
        configd: dict[str, Any],
        tagger_pool: list[type[ReadAwareProcess]],
        flagger_pool: list[type[ReadAwareProcess]]
    ) -> Self:
        if not set(['read-processors', 'variant-flaggers', 'opts']) == configd.keys():
            raise ValueError('config does not have correct top-level keys')

        excludes_opt= cast(bool, configd['opts']['mandate-excludes'])
        tagger_paramd = cast(dict[str, Any], configd['read-processors'])
        flagger_paramd = cast(dict[str, Any], configd['variant-flaggers'])
        if not isinstance(excludes_opt, bool):
            raise ValueError
        if not isinstance(tagger_paramd, dict):
            raise ValueError
        if not isinstance(flagger_paramd, dict):
            raise ValueError

        ra_ex = cls(None, excludes_opt)
        if not ra_ex.check_namespacing_clashfree(tagger_pool + flagger_pool):
            raise ValueError('Processes with clashing namespaces')

        tagger_nsd = {proc_type.ProcessNamespace: proc_type for proc_type in tagger_pool}
        flagger_nsd = {proc_type.ProcessNamespace: proc_type for proc_type in flagger_pool}

        if not set(tagger_paramd.keys()).issubset(tagger_nsd):
            raise ValueError('config requesting unknown taggers')
        if not set(flagger_paramd.keys()).issubset(flagger_nsd):
            raise ValueError('config requesting unknown flaggers')

        # TODO: if using pydantic doesn't engender end user having to necessarily define paramater dataclasses,
        # then clearly need pydantic validated config dataclass (may suffice in the meantime anyway since this is hp2 not htsflow)
        # TODO: stop using hypens, given config keys will mostly correspond to python symbols
        # config order defines execution order
        ordered_req_taggers = [tagger_nsd[proc_name](config_params, config_params['require-marks'], config_params['exclude-marks']) for proc_name, config_params in tagger_paramd.items()]
        ordered_req_flaggers = [flagger_nsd[proc_name](config_params, config_params['require-marks'], config_params['exclude-marks']) for proc_name, config_params in flagger_paramd.items()]

        _, tags_added = ra_ex.validate_tagger_order(ordered_req_taggers, ra_ex.mandate_excludes, raise_on_fail=True)
        _ = ra_ex.validate_flagger_exec(ordered_req_flaggers, tags_added, ra_ex.mandate_excludes, raise_on_fail=True)

        ra_ex._taggers = tuple(ordered_req_taggers)
        ra_ex._flaggers = tuple(ordered_req_flaggers)

        return ra_ex

    def run(  # TODO
        self,
        run_data: RunParams
    ):
        for proc in self.taggers:
            proc(run_data)  # TODO: can now fix typing of these classes
            # NOTE: remember decorators are busted

    

