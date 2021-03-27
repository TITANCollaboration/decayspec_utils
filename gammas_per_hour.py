import pandas as pd

class gammas:
    def __init__(self, filename, output_filename, initial_decays, header):
        self.filename = filename
        self.initial_decays = initial_decays  # for one hour sb129m1
        #self.start_header = "sb129_1850_bg"
        self.start_header = header #"sb129m1_bg"
        self.output_filename = output_filename
        #self.initial_decays = 27  # for one hour sb129
        #self.start_header = "sb129_bg"

        return

    def write_geant4_macro(self, gamma_energy, probability_of_decay):
        with open(self.output_filename, 'a') as geantMacroFile:
            num_decays = round((probability_of_decay / 100) * (self.initial_decays))
            if num_decays > 0:
                geantMacroFile.write("# START: %s_%s\n" % (self.start_header, gamma_energy))
                geantMacroFile.write("/pga/selectGunAction 1\n")
                geantMacroFile.write("/gps/particle gamma\n")
                geantMacroFile.write("/gps/position 0 0 0\n")
                geantMacroFile.write("/gps/energy %s keV\n" % (gamma_energy))
                geantMacroFile.write("/gps/ang/type iso\n")
                geantMacroFile.write("/histo/filename %s_%s\n" % (self.start_header, gamma_energy))
                geantMacroFile.write("/run/beamOn %s\n" % (num_decays))
                geantMacroFile.write("# END\n")

# END

    def read_in_gammas(self):
        df = pd.read_csv(self.filename, sep=',')
        #print(df)
        for i, gamma_info in df.iterrows():
            #print(gamma_info[0], gamma_info[1])
            #print(gamma_info[0])
            self.write_geant4_macro(gamma_info[0], gamma_info[1])
        return


#mygamma = gammas("gammas_sb129m1.csv", "sb129_g4_background.mac", 460847, "sb129m1_bg")
mygamma = gammas("gammas_sb129.csv", "sb129_g4_background.mac", 1155, "sb129_bg")

mygamma.read_in_gammas()
