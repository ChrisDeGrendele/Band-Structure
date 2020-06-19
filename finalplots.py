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
save_dir = "newdata/png/"
dpi_rez = 300



def doplots(file_list):
    for file in file_list:

        data = scp.loadmat(data_dir + file) #get data
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

        plt.savefig(save_dir + "{}.png".format(file[:-4]), dpi=dpi_rez)
        print(file)

def do_individual_vertical_bands(file_list):
    #plt.figure()
    i = 1
    for file in file_list:

        data = scp.loadmat(data_dir + file) #get data
        Imag_E = data['Imag_E'].flatten()
        Real_E = data['Real_E'].flatten()
        Imag_K = data['Imag_k'].flatten()
        Real_K = data['Real_k'].flatten()

        plt.subplot(220+i)
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

        i+=1

        plt.subplots_adjust(wspace=0,hspace=0)
        if not os.path.exists('png'):
            os.makedirs('png')

        plt.savefig(save_dir + "{}.png".format(file[:-4]), dpi=dpi_rez)
        print(file)

def do_vertical_bands(file_list):
    i = 0
    fig = plt.figure(figsize=(10, 8))
    outer = gridspec.GridSpec(2, 2, wspace=0.2, hspace=0.2)

    for file in file_list:

        data = scp.loadmat(data_dir + file) #get data
        Imag_E = data['Imag_E'].flatten()
        Real_E = data['Real_E'].flatten()
        Imag_K = data['Imag_k'].flatten()
        Real_K = data['Real_k'].flatten()

        inner = gridspec.GridSpecFromSubplotSpec(1,2,
                    subplot_spec=outer[i],wspace=0,hspace=0)

        #Imaginary Data:
        # I think the scatter plot looks nicer on imaginary side
        #ax1.plot(np.abs(Imag_K), Imag_E,'.', color="b", markersize=1)

        ax1 = plt.Subplot(fig, inner[0])
        ax1.plot(np.abs(Imag_K), Imag_E,'.', color="b", markersize=0.5)
        ax1.set_xlim(np.nanmax(Imag_K[Imag_K != np.inf])+0.5,0)
        ax1.set_ylabel(r'$Energy \enspace (eV)$')#, fontsize=12)
        ax1.set_xlabel(r'$Imaginary \enspace \vec{k}$')
        ax1.set_ylim((-5,5))
        fig.add_subplot(ax1)

        #Real Data:
        #Cut data into pos/neg and connect lines that way to avoid really weird zigzag plots

        Pos_Real_K = Real_K[Real_E > 0]
        Pos_Real_E = Real_E[Real_E > 0]
        Neg_Real_K = Real_K[Real_E < 0]
        Neg_Real_E = Real_E[Real_E < 0]

        ax2 = plt.Subplot(fig, inner[1])
        ax2.plot(np.abs(Pos_Real_K), Pos_Real_E,'.', color="b", markersize=0.5)
        ax2.plot(np.abs(Neg_Real_K), Neg_Real_E,'.', color="b", markersize=0.5)
        ax2.set_xticks([1,2,3])
        ax2.set_yticks([])
        ax2.set_ylim((-5,5))
        ax2.set_xlim(0,np.nanmax(Real_K[Real_K != np.inf]))
        ax2.set_xlabel(r'$Real \enspace \vec{k}$')
        fig.add_subplot(ax2)
        #f.text(0.06, 0.5, r'$Energy \enspace (eV)$', ha='center', va='center', rotation='vertical')
        #f.text(0.333, 0.03, r'$Imaginary \enspace \vec{k}$', ha='center', va='center', rotation='horizontal')
        #f.text(0.7, 0.03, r'$Real \enspace \vec{k}$', ha='center', va='center', rotation='horizontal')
        #f.text(0.5, 0.94, file[:-4], ha='center', va='center', rotation='horizontal')

        i+=1

        #plt.subplots_adjust(wspace=0,hspace=0)
        #if not os.path.exists('png'):
        #    os.makedirs('png')
        #plt.savefig("png/{}.png".format(file[:-4]), dpi=dpi_rez)
        ##print(file)
    plt.savefig(save_dir + "combined_vertical_bands.png", dpi=dpi_rez)

def do_benzene(file_list):
    #i = 0

    fig = plt.figure(figsize=(10, 8))
    outer = gridspec.GridSpec(13,13, wspace=0.2, hspace=0.2)

    for file in file_list:

        data = scp.loadmat(data_dir + file) #get data
        Imag_E = data['Imag_E'].flatten()
        Real_E = data['Real_E'].flatten()
        Imag_K = data['Imag_k'].flatten()
        Real_K = data['Real_k'].flatten()


        if file == 'Meta.mat':
            inner = gridspec.GridSpecFromSubplotSpec(1,2,
                    subplot_spec=outer[2:6, 2:6],wspace=0,hspace=0)
        elif file == 'Ortho.mat':
            inner = gridspec.GridSpecFromSubplotSpec(1,2,
                    subplot_spec=outer[2:6, 8:12],wspace=0,hspace=0)
        elif file == 'Para.mat':
            inner = gridspec.GridSpecFromSubplotSpec(1,2,
                    subplot_spec=outer[7:11,5:9],wspace=0,hspace=0)



        #Imaginary Data:
        # I think the scatter plot looks nicer on imaginary side
        #ax1.plot(np.abs(Imag_K), Imag_E,'.', color="b", markersize=1)

        ax1 = plt.Subplot(fig, inner[0])
        ax1.plot(np.abs(Imag_K), Imag_E,'.', color="b", markersize=0.5)
        ax1.set_xlim(np.nanmax(Imag_K[Imag_K != np.inf])+0.5,0)
        ax1.set_ylabel(r'$Energy \enspace (eV)$')#, fontsize=12)
        ax1.set_xlabel(r'$Imaginary \enspace \vec{k}$')
        if file == "Para.mat":
            ax1.set_xticks([0,2.5,5])
        else:
            ax1.set_xticks([0,10,20,30])

        ax1.set_ylim((-5,5))


        fig.add_subplot(ax1)




        #Real Data:
        #Cut data into pos/neg and connect lines that way to avoid really weird zigzag plots

        Pos_Real_K = Real_K[Real_E > 0]
        Pos_Real_E = Real_E[Real_E > 0]
        Neg_Real_K = Real_K[Real_E < 0]
        Neg_Real_E = Real_E[Real_E < 0]

        ax2 = plt.Subplot(fig, inner[1])
        ax2.plot(np.abs(Pos_Real_K), Pos_Real_E, color="b", markersize=0.5)
        ax2.plot(np.abs(Neg_Real_K), Neg_Real_E, color="b", markersize=0.5)
        ax2.set_xticks([1,2,3])
        ax2.set_yticks([])
        ax2.set_ylim((-5,5))
        ax2.set_xlim(0,np.nanmax(Real_K[Real_K != np.inf]))
        ax2.set_xlabel(r'$Real \enspace \vec{k}$', fontsize=12)
        fig.add_subplot(ax2)
        #i+=2

    plt.savefig(save_dir + "combined_benzene.png", dpi=dpi_rez, bbox_inches = 'tight',pad_inches = 0.2)

file_list = ['simple2bandex.mat', 'CombMolecule.mat', 'Meta.mat',
            'Ortho.mat', 'Para.mat', 'siloxane2.mat']#, 'Azulene.mat', 'Assymetric.mat']
vertical_bands = ['VerticalBand.1.mat', 'VerticalBand.01.mat','VerticalBand.001.mat',
                  'VerticalBand.0001.mat']
benzene = ['Meta.mat','Ortho.mat', 'Para.mat']


#New files generated by new matlab script... For comparing results, to see if it
#made the same figures.



doplots(file_list)
do_vertical_bands(vertical_bands)
#do_individual_vertical_bands(vertical_bands)
do_benzene(benzene)
