import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from datetime import datetime
start_time = datetime.now()

def analysis(typ, endresidue): #specify type of glycosilation (in string) and end residue# until to analyze
    fig1 = plt.figure() #rmsd
    ax1 = fig1.add_subplot(111)

    fig2 = plt.figure() #rgyr
    ax2 = fig2.add_subplot(111)

    traj, protein_sel, protein_traj, rmsd, rmsdfixed, rgyr, rgyrfixed = [None] * 6, [None] * 6, [None] * 6, [None] * 6, [None] * 6, [None] * 6, [None] * 6
    traj[1] = md.load('sosip_'+typ+'/amber_run/01'+typ+'.nc', top='sosip_'+typ+'/amber_run/01'+typ+'.pdb')
    protein_sel[1] = traj[1].topology.select('protein and backbone')
    protein_traj[1] = traj[1].atom_slice(protein_sel[1])
    rmsd[1] = md.rmsd(protein_traj[1], protein_traj[1], 0)
    rmsdfixed[1] = rmsd[1][1:endresidue]
    ax1.plot(np.arange(0, len(rmsdfixed[1])), rmsdfixed[1]*10, color='b', label='01'+typ)
    rgyr[1] = md.compute_rg(protein_traj[1])
    rgyrfixed[1] = rgyr[1][1:endresidue]
    ax2.plot(np.arange(0, len(rgyrfixed[1])), rgyrfixed[1]*10, color='b', label='01'+typ)
    colors = "grcmykw"
    color_index = 0
    for i in range(2,6):
        traj[i] = md.load('0'+str(i)+typ+'/amber_run/0'+str(i)+typ+'.nc', top='0'+str(i)+typ+'/amber_run/0'+str(i)+typ+'.pdb')
        protein_sel[i] = traj[i].topology.select('protein and backbone')
        protein_traj[i] = traj[i].atom_slice(protein_sel[i])
        rmsd[i] = md.rmsd(protein_traj[i], protein_traj[i], 0)
        rmsdfixed[i] = rmsd[i][1:endresidue]
        ax1.plot(np.arange(0, len(rmsdfixed[i])), rmsdfixed[i]*10, c=colors[color_index], label='0'+str(i)+typ)
        rgyr[i] = md.compute_rg(protein_traj[i])
        rgyrfixed[i] = rgyr[i][1:endresidue]
        ax2.plot(np.arange(0, len(rgyrfixed[i])), rgyrfixed[i]*10, c=colors[color_index], label='0'+str(i)+typ)
        color_index += 1

     #print(rmsd1fixed)
    ax1.set_xlabel('Time (ns)')
    ax1.set_ylabel('RMSD (A)')
    ax1.legend()
    fig1.savefig('Analysis/rmsd_'+typ+'.png', bbox_inches='tight', pad_inches=0.1, dpi = 480, facecolor='white' )
    print('Analysis/rmsd_'+typ+'.png file generated! Last updated at ',datetime.now())
    ax2.set_xlabel('Time (ns)')
    ax2.set_ylabel('rgyr (A)')
    ax2.legend()
    fig2.savefig('Analysis/rgyr_'+typ+'.png', bbox_inches='tight', pad_inches=0.1, dpi = 480, facecolor='white' )
    print('Analysis/rgyr_'+typ+'.png file generated! Last updated at ',datetime.now())


analysis('prot', 500)
analysis('m9', 300)
analysis('nat', 200)
end_time = datetime.now()
print('Program took {}'.format(end_time - start_time))
