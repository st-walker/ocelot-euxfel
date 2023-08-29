from collections import defaultdict


class PlacedPhysicsProcess:
    def __init__(self, proc, start, stop, description="", active=True):
        self.proc = proc
        self.start_name = start
        self.stop_name = stop
        self.active = active


class PlacedElementProcess(PlacedPhysicsProcess):
    def __init__(self, proc, start, stop, description="", active=True):
        self.proc = proc
        self.start_name = start
        self.stop_name = stop
        self.active = active

    def modify_element(self, element):
        pass

    def modify_sequence(self, element):
        pass


class LaserHeaterProcess(PlacedElementProcess):
    def __init__(self, proc, start, stop, undulator_name, description=""):
        pass


class PhysicsList:
    def __init__(self, processes):
        self._processes = defaultdict(list)
        for process in processes:
            processes[(type(process).__name__)].append(process)

    def add_to_navi(self, navi):
        pass

    def __getitem__(self, phys):
        try:
            return self._processes[phys]
        except KeyError:
            pass

        try:
            return self._processes[phys.__name__]
        except (KeyError, AttributeError):
            pass

        raise KeyError("Physics process not found")

    def add_physics_process(self, process, start, stop):
        pass
