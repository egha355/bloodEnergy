class Problem_Params:
    def __init__(self):
        self.timeIncrement = 0.01
        self.startTime = 0.0
        self.timeSteps = 301 #Hi
        self.outputFrequency = 1
        # properties   
        self.diffusivity = 1.57e-7
     #   self.conductivity = 0.42
     #   self.rhoC = 4.0e6 
     #   self.source = 0.0
     #   self.convection = 20
     #   self.Tair = 0.0
     #   self.velocity = 0.104562
        self.Tinit = 37.0
        self.elementsFile = "input/elements.csv"
        self.nodesFile = "input/nodes.csv"
        self.Tinlet = 37.0

    def time_increment_set(self, value):
        self.timeIncrement = value
