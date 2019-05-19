import pint

class GlobalRegistry():
    __instance = None
    @staticmethod 
    def getRegistry():
        """ Static access method. """
        if GlobalRegistry.__instance == None:
            GlobalRegistry()
        return GlobalRegistry.__instance

    def __init__(self):
        """ Virtually private constructor. """
        if GlobalRegistry.__instance != None:
            raise Exception("This class is a singleton!")
        else:
            self.ureg = pint.UnitRegistry()
            GlobalRegistry.__instance = self