import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.io as scp
import numpy as np
import os

file_list = ['simple2bandex.mat', 'CombMolecule.mat', 'Meta.mat',
            'Ortho.mat', 'Para.mat', 'siloxane2.mat', 'Assymetric.mat']#, 'Azulene.mat']
vertical_bands = ['VerticalBand.1.mat', 'VerticalBand.01.mat','VerticalBand.001.mat',
                  'VerticalBand.0001.mat']
benzene = ['Meta.mat','Ortho.mat', 'Para.mat']

data_dir = "newdata/"
save_dir = "newdata/png"
dpi_rez = 300



def doplots(file_list):
    for file in file_list:

        data = scp.loadmat(new_dir + file) #get data
        Imag_E = data['Imag_E'].flatten()
        Real_E = data['Real_E'].flatten()
        Imag_K = data['Imag_k'].flatten()
        Real_K = data['Real_k'].flatten()

        f, (ax1,ax2) = plt.subplots(1,2,sharey=True)

        #Imaginary Data:
        # I think the scatter plot looks nicer on imaginary side
        ax1.plot(np.abs(Imag_K), Imag_E,'.', color="b", markersize=1)

        #Real Data:

        #Cut data into pos/neg and connect lines that way to avoid really weird zigzag plots
        if (file == "siloxane2.mat") or (file == "silo_highrez.mat") or (file=="Assymetric.mat") : #problem system
            ax2.plot(np.abs(Real_K), Real_E,'.', color="b", markersize=1)
        else:
            Pos_Real_K = Real_K[Real_E > 0]
            Pos_Real_E = Real_E[Real_E > 0]
            Neg_Real_K = Real_K[Real_E < 0]
            Neg_Real_E = Real_E[Real_E < 0]


            ax2.plot(np.abs(Pos_Real_K), Pos_Real_E,'-', color="b", markersize=1)
            ax2.plot(np.abs(Neg_Real_K), Neg_Real_E,'-', color="b", markersize=1)

        ax1.set_xlim(np.nanmax(Imag_K[Imag_K != np.inf])+0.5,0)

        ax2.set_xlim(0,np.nanmax(Real_K[Real_K != np.inf]))
        f.text(0.06, 0.5, r'$Energy \enspace (eV)$', ha='center', va='center', rotation='vertical')
        f.text(0.333, 0.03, r'$Imaginary \enspace \vec{k}$', ha='center', va='center', rotation='horizontal')
        f.text(0.7, 0.03, r'$Real \enspace \vec{k}$', ha='center', va='center', rotation='horizontal')
        #f.text(0.5, 0.94, file[:-4], ha='center', va='center', rotation='horizontal')

        plt.subplots_adjust(wspace=0,hspace=0)
        if not os.path.exists('png'):
            os.makedirs('png')

        if new_files:
            save_dir = "../newdata/png/"
        else:
            save_dir = "png/"

        plt.savefig(save_dir + "{}.png".format(file[:-4]), dpi=dpi_rez)
        print(file)
        plt.show()
