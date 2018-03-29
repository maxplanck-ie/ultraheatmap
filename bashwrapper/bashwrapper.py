import subprocess

pipeout = []

class Bash():
    def __init__(self, command):
        self.command = command
        self.__stderr = None
        self.__stdout = None
        self.p_status = None

        self.proc = subprocess.Popen(self.command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        self.__stdout, self.__stderr = self.proc.communicate()
        self.p_status = self.proc.wait()
        if self.p_status != 0:

            raise Exception(self.__stderr)
            self.proc.kill()

    def stderr(self):
        return(self.__stderr)

    def stdout(self):
        return(self.__stdout)
