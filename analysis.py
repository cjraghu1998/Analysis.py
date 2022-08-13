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
    
    fig3 = plt.figure() #rmsf 1-633
    ax3 = fig3.add_subplot(111)
    
    fig4 = plt.figure() #rmsf 634:1267
    ax4 = fig4.add_subplot(111)
    
    fig5 = plt.figure() #rmsf 1268:1901
    ax5 = fig5.add_subplot(111)

    rmsdrange, rgyrrange, rmsfrange, rmsfrange2, rmsfrange3 = [None] * 5, [None] * 5, [None] * 5, [None] * 5, [None] * 5
    traj = md.load('sosip_prot/amber_run/01prot.nc', top='sosip_prot/amber_run/01prot.pdb') #0th frame of prot
    protein_sel = traj.topology.select('protein and name CA')
    protein_traj0 = traj.atom_slice(protein_sel)
    
    colors = "grcmykw"
    color_index = 0
    for i in range(1,6):
        traj = md.load('0'+str(i)+typ+'/amber_run/0'+str(i)+typ+'.nc', top='0'+str(i)+typ+'/amber_run/0'+str(i)+typ+'.pdb')
        protein_sel = traj.topology.select('protein and name CA')
        protein_traj = traj.atom_slice(protein_sel)
        rmsd = md.rmsd(protein_traj, protein_traj, 0)
        rmsdrange[i-1] = rmsd[1:endresidue]
        ax1.plot(np.arange(0, len(rmsdrange[i-1])), rmsdrange[i-1]*10, c=colors[color_index], label='0'+str(i)+typ)
        rgyr = md.compute_rg(protein_traj)
        rgyrrange[i-1] = rgyr[1:endresidue]
        ax2.plot(np.arange(0, len(rgyrrange[i-1])), rgyrrange[i-1]*10, c=colors[color_index], label='0'+str(i)+typ)
        rmsf = md.rmsf(protein_traj, protein_traj0, 0)
        rmsfrange[i-1] = rmsf[0:633]
        rmsfrange2[i-1] = rmsf[634:1267]
        rmsfrange3[i-1] = rmsf[1268:1901]
        color_index += 1


    rmsfrange10 = np.multiply(rmsfrange,10)
    rmsfrange20 = np.multiply(rmsfrange2,10)
    rmsfrange30 = np.multiply(rmsfrange3,10)

    ax3.plot(np.arange(0, len(rmsfrange[0])), np.mean(rmsfrange10, axis=0), color='b', label="mean", linewidth=0.7)

    ax3.fill_between(np.arange(0, len(rmsfrange[0])), np.mean(rmsfrange10, axis=0) - np.std(rmsfrange10, axis=0), np.mean(rmsfrange10, axis=0) + np.std(rmsfrange10, axis=0), color='b', alpha=0.4)
    
    ax4.plot(np.arange(0, len(rmsfrange2[0])), np.mean(rmsfrange20, axis=0), color='b', label="mean", linewidth=0.7)

    ax4.fill_between(np.arange(0, len(rmsfrange2[0])), np.mean(rmsfrange20, axis=0) - np.std(rmsfrange20, axis=0), np.mean(rmsfrange20, axis=0) + np.std(rmsfrange20, axis=0), color='b', alpha=0.4)
    
    ax5.plot(np.arange(0, len(rmsfrange3[0])), np.mean(rmsfrange30, axis=0), color='b', label="mean", linewidth=0.7)

    ax5.fill_between(np.arange(0, len(rmsfrange3[0])), np.mean(rmsfrange30, axis=0) - np.std(rmsfrange30, axis=0), np.mean(rmsfrange30, axis=0) + np.std(rmsfrange30, axis=0), color='b', alpha=0.4)

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
    ax3.set_xlabel('Residues')
    ax3.set_ylabel('rmsf (A)')
    ax3.legend()
    fig3.savefig('Analysis/rmsf1_'+typ+'.png', bbox_inches='tight', pad_inches=0.1, dpi = 480, facecolor='white' )
    print('Analysis/rmsf1_'+typ+'.png file generated! Last updated at ',datetime.now())
    ax4.set_xlabel('Residues')
    ax4.set_ylabel('rmsf (A)')
    ax4.legend()
    fig4.savefig('Analysis/rmsf2_'+typ+'.png', bbox_inches='tight', pad_inches=0.1, dpi = 480, facecolor='white' )
    print('Analysis/rmsf2_'+typ+'.png file generated! Last updated at ',datetime.now())
    ax5.set_xlabel('Residues')
    ax5.set_ylabel('rmsf (A)')
    ax5.legend()
    fig5.savefig('Analysis/rmsf3_'+typ+'.png', bbox_inches='tight', pad_inches=0.1, dpi = 480, facecolor='white' )
    print('Analysis/rmsf3_'+typ+'.png file generated! Last updated at ',datetime.now())


analysis('prot', 500)
analysis('m9', 300)
analysis('nat', 150)
end_time = datetime.now()
print('Program took {}'.format(end_time - start_time))
