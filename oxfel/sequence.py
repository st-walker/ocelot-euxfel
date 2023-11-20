from __future__ import annotations

from collections.abc import Sequence
from collections import Counter
from typing import Optional, Union, Type, Iterable, TypeVar, Callable, Any

from ocelot.cpbd.magnetic_lattice import flatten

from ocelot.cpbd.elements.optic_element import OpticElement

ElementT = TypeVar("ElementT", bound=OpticElement)
ElementAccessType = Optional[Union[int, str, ElementT]]


class ElementAccessError(KeyError):
    pass

class MachineSequence(Sequence):
    def __init__(self, sequence: list[ElementT], *args, **kwargs):
        """Name addressable sequence of elements"""
        super().__init__(*args, **kwargs)
        self._sequence = list(flatten(sequence))
        name_count = Counter(self.names())
        if not set(name_count.values()).issubset({1}):
            raise ValueError("Duplicate names found")

    def __getitem__(self, key: ElementAccessType) -> Union[ElementT, MachineSequence]:
        if isinstance(key, slice):
            start, step, stop = key.start, key.step, key.stop
            if step is not None:
                raise TypeError("Slice step is not supported")

            if start is None:
                start = 0

            start = self._normalise_key(start)
            stop = self._normalise_key(stop)
            return MachineSequence(self._sequence[start:stop])
        else:
            index = self._normalise_key(key)
            return self._sequence[index]

    def closed_interval(
        self, start: ElementAccessType = None, stop: ElementAccessType = None
    ) -> MachineSequence:
        """Get interval including the stop element at the end (unlike __getitem__)"""
        if start is None:
            start = 0
        istart = self._normalise_key(start)
        if stop is None:
            stop = -1
        istop = self._normalise_key(stop)
        # Closed interval so add 1:
        istop += 1
        return self[istart:istop]

    def __iter__(self):
        yield from iter(self._sequence)

    def _normalise_key(self, key: ElementAccessType) -> int:
        # If already an int then just return it
        if key is None:
            return None

        if isinstance(key, int):
            if key >= 0:
                return key
            # Convert negative index into a positive one.
            return key % len(self)

        # If a name of an element find the index
        if isinstance(key, str):
            return self.names().index(key)

        # if an element instance then get the key
        try:
            return self._sequence.index(key)
        except (TypeError, ValueError):
            pass

        raise ElementAccessError(f"Unable to normalise key: {key}")

    def __len__(self) -> int:
        return len(self._sequence)

    def __add__(self, other: Iterable[ElementT]) -> MachineSequence:
        return type(self)(list(self) + list(flatten(other)))

    def __str__(self) -> str:
        strs = "\n".join([repr(s) for s in self._sequence])
        return f"<{type(self).__name__}:\n{strs}>"

    def __contains__(self, key: ElementAccessType) -> bool:
        try:
            return key in self._sequence or key.id in self.names()
        except AttributeError:
            return key in self.names()

    def names(self) -> list[str]:
        return [x.id for x in self]

    def element_attributes(self, key, property_name: str) -> np.array:
        indices = [i for (i, ele) in enumerate(self) if (ele.id == key) or (ele is key)]
        if not indices:
            raise KeyError
        result = []
        for i in indices:
            result.append(getattr(self[i], property_name))
        return np.squeeze(np.array(result))

    def length(self) -> float:
        return sum(x.l for x in self)

    def element_s(self, key: Union[str, ElementT]) -> float:
        indices = [i for (i, ele) in enumerate(self) if (ele.id == key) or (ele is key)]

        if not indices:
            raise KeyError

        result = []
        s = 0
        for i, element in enumerate(self):
            length = element.l
            s += length
            if i in indices:
                result.append(s)
        return np.squeeze(np.array(result))

    def __repr__(self) -> str:
        return f"<MachineSequence: {repr(self._sequence)}>"

    def reverse(self):
        """reverse *IN PLACE*"""
        self._sequence.reverse()
