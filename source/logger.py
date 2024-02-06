import os


class Logger:
    def __init__(self, path=".", name="log.txt"):
        if not name.endswith(".txt"):
            name += ".txt"
        self.logpath = os.path.join(path, name)
        self._empty = True

    def write(self, line):
        if self._empty:
            with open(self.logpath, "w+") as self.logfile:
                self._write(line)
            self._empty = False
        else:
            with open(self.logpath, "a") as self.logfile:
                self._write(line)

    def _write(self, line):
        if isinstance(line, list):
            for w in line:
                self.logfile.write(str(w) + "\n")
        else:
            self.logfile.write(str(line) + "\n")
