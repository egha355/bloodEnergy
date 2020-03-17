class Problem_Params:
    def __init__(self):
        self.timeIncrement = 10.0
        self.startTime = 0.0
        self.timeSteps = 101 #Hi
        self.diffusionOutputFrequency = 100
        self.conductivity = 0.42
        self.rhoC = 4.0e6 
        self.source = 0.0
        self.convection = 20
        self.Tair = 0.0
        self.Tinit = 37.0
        self.tissueElementsFile = "input/cylinder_elements.csv"
        self.tissueNodesFile = "input/cylinder_nodes.csv"
        self.numberOfNodes = 11476
        self.numberOfElements = 49243


    def time_increment_set(self, value):
        self.timeIncrement = value
