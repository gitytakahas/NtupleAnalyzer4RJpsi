import os, math, sys
from ROOT import TFile, TH1F, gROOT, TTree, Double, TChain
import numpy as num

gROOT.SetBatch(True)

from optparse import OptionParser, OptionValueError
usage = "usage: python runTauDisplay.py"
parser = OptionParser(usage)

parser.add_option("-o", "--out", default='Myroot.root', type="string", help="output filename", dest="out")
parser.add_option("-p", "--path", default='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi_20191002_BcJpsiMuNu_020519/BcJpsiMuNu_020519/BcJpsiMuNu_020519_v1/191002_132739/0000/', type="string", help="path", dest="path")

(options, args) = parser.parse_args()

print 'output file name = ', options.out

outputfile = TFile('root/' + options.out, 'recreate')
otree = TTree('tree', 'tree')

#file2include = 'dcap://t3se01.psi.ch:22125/' + options.path + '/*.root'
#file2include = 'dcap://t3se01.psi.ch:22125/' + options.path + '/*.root'
#file2include = 'root://storage01.lcg.cscs.ch/pnfs/lcg.cscs.ch/cms/trivcat/store/user/ytakahas//RJpsi_20191031_v3_BcJpsiTauNu_020519/BcJpsiTauNu_020519/BcJpsiTauNu_020519_v1/191031_113035/0000/*.root'
file2include = 'root://storage01.lcg.cscs.ch/pnfs/lcg.cscs.ch/cms/trivcat/store/user/ytakahas/'
#RJpsi_2019-11-01-160633_BcJpsiTauNu_020519/BcJpsiTauNu_020519/BcJpsiTauNu_020519dz_zero/191101_150803/0000/*.root'

#file2include += 'RJpsi_2019-11-04-094448_BcJpsiTauNu_020519/BcJpsiTauNu_020519/BcJpsiTauNu_020519_20191104094448_dz0/191104_084614/0000/*.root'
#file2include += 'RJpsi_2019-11-04-112306_BcJpsiTauNu_020519/BcJpsiTauNu_020519/BcJpsiTauNu_020519_20191104112306_pv0/191104_102445/0000/*.root'
#file2include += 'RJpsi_2019-11-04-111601_BcJpsiTauNu_020519//BcJpsiTauNu_020519/BcJpsiTauNu_020519_20191104111601_dz0p5/191104_101731/0000/*.root'
#file2include += 'RJpsi_2019-11-04-111918_BcJpsiTauNu_020519/BcJpsiTauNu_020519/BcJpsiTauNu_020519_20191104111918_pv0p5/191104_102051/0000/*.root'
file2include += options.path
file2include += '/*.root'



print 'file2include = ', file2include 

chain = TChain('ntuplizer/tree', 'tree')

chain.Add(file2include)

chain.SetBranchStatus('*', 0)
chain.SetBranchStatus('EVENT_*', 1)
chain.SetBranchStatus('JpsiTau_tau*', 1)
chain.SetBranchStatus('JpsiTau_B*', 1)
chain.SetBranchStatus('JpsiTau_bbPV*', 1)
chain.SetBranchStatus('JpsiTau_isgen3matched', 1)
chain.SetBranchStatus('JpsiTau_isgen3', 1)
chain.SetBranchStatus('JpsiTau_ngentau', 1)
chain.SetBranchStatus('JpsiTau_isRight*', 1)
chain.SetBranchStatus('JpsiTau_pf*', 1)

#otree = chain.CloneTree(0)


#B_mass = num.zeros(1, dtype=float)
#B_pt = num.zeros(1, dtype=float)
#otree.Branch('B_mass', B_mass, 'B_mass/D')
#otree.Branch('B_pt', B_pt, 'B_pt/D')


tau_pt = num.zeros(1, dtype=float)
tau_eta = num.zeros(1, dtype=float)
tau_phi = num.zeros(1, dtype=float)
tau_mass = num.zeros(1, dtype=float)
tau_vx = num.zeros(1, dtype=float)
tau_vy = num.zeros(1, dtype=float)
tau_vz = num.zeros(1, dtype=float)

tau_lip = num.zeros(1, dtype=float)
tau_lips = num.zeros(1, dtype=float)
tau_pvip = num.zeros(1, dtype=float)
tau_pvips = num.zeros(1, dtype=float)
tau_fl3d = num.zeros(1, dtype=float)
tau_fls3d = num.zeros(1, dtype=float)
tau_alpha = num.zeros(1, dtype=float)

tau_max_dr_3prong = num.zeros(1, dtype=float)
#tau_iso03 = num.zeros(1, dtype=float)
tau_charge = num.zeros(1, dtype=int)
tau_vprob = num.zeros(1, dtype=float)
#tau_dz2pv = num.zeros(1, dtype=float)
#tau_dxy2pv = num.zeros(1, dtype=float)
tau_isRight = num.zeros(1, dtype=bool)
tau_isgen3matched = num.zeros(1, dtype=int)
tau_isgen3 = num.zeros(1, dtype=int)
#tau_ngentau = num.zeros(1, dtype=int)
#tau_nch = num.zeros(1, dtype=int)
tau_gentaupt = num.zeros(1, dtype=float)
#
b_mass = num.zeros(1, dtype=float)
b_pt = num.zeros(1, dtype=float)
b_eta = num.zeros(1, dtype=float)
b_phi = num.zeros(1, dtype=float)

tau_isRight_bS = num.zeros(1, dtype=bool)
tau_isRight_aS = num.zeros(1, dtype=bool)
tau_isRight_bS_ith = num.zeros(1, dtype=int)
tau_isRight_aS_ith = num.zeros(1, dtype=int)
tau_isRight_bS_n = num.zeros(1, dtype=int)
tau_isRight_aS_n = num.zeros(1, dtype=int)

event = num.zeros(1, dtype=int)
run = num.zeros(1, dtype=int)
lumi = num.zeros(1, dtype=int)

vx = num.zeros(1, dtype=float)
vy = num.zeros(1, dtype=float)
vz = num.zeros(1, dtype=float)
#
#jpsi_pt = num.zeros(1, dtype=float)
#jpsi_eta = num.zeros(1, dtype=float)
#jpsi_phi = num.zeros(1, dtype=float)
#jpsi_mass = num.zeros(1, dtype=float)
#jpsi_vprob = num.zeros(1, dtype=float)
#jpsi_flightSig3D = num.zeros(1, dtype=float)
#
#
otree.Branch('tau_pt', tau_pt, 'tau_pt/D')
otree.Branch('tau_eta', tau_eta, 'tau_eta/D')
otree.Branch('tau_phi', tau_phi, 'tau_phi/D')
otree.Branch('tau_mass', tau_mass, 'tau_mass/D')
otree.Branch('tau_vx', tau_vx, 'tau_vx/D')
otree.Branch('tau_vy', tau_vy, 'tau_vy/D')
otree.Branch('tau_vz', tau_vz, 'tau_vz/D')

otree.Branch('tau_lip', tau_lip, 'tau_lip/D')
otree.Branch('tau_lips', tau_lips, 'tau_lips/D')
otree.Branch('tau_pvip', tau_pvip, 'tau_pvip/D')
otree.Branch('tau_pvips', tau_pvips, 'tau_pvips/D')
otree.Branch('tau_fl3d', tau_fl3d, 'tau_fl3d/D')
otree.Branch('tau_fls3d', tau_fls3d, 'tau_fls3d/D')
otree.Branch('tau_alpha', tau_alpha, 'tau_alpha/D')


otree.Branch('tau_max_dr_3prong', tau_max_dr_3prong, 'tau_max_dr_3prong/D')
#otree.Branch('tau_iso03', tau_iso03, 'tau_iso03/D')
otree.Branch('tau_charge', tau_charge, 'tau_charge/I')
otree.Branch('tau_vprob', tau_vprob, 'tau_vprob/D')
#otree.Branch('tau_dz2pv', tau_dz2pv, 'tau_dz2pv/D')
#otree.Branch('tau_dxy2pv', tau_dxy2pv, 'tau_dxy2pv/D')
otree.Branch('tau_isRight', tau_isRight, 'tau_isRight/O')
otree.Branch('tau_isgen3matched', tau_isgen3matched, 'tau_isgen3matched/I')
otree.Branch('tau_isgen3', tau_isgen3, 'tau_isgen3/I')
#otree.Branch('tau_ngentau', tau_ngentau, 'tau_ngentau/I')
#otree.Branch('tau_nch', tau_nch, 'tau_nch/I')
otree.Branch('tau_gentaupt', tau_gentaupt, 'tau_gentaupt/D')
#
otree.Branch('b_mass', b_mass, 'b_mass/D')
otree.Branch('b_pt', b_pt, 'b_pt/D')
otree.Branch('b_eta', b_eta, 'b_eta/D')
otree.Branch('b_phi', b_phi, 'b_phi/D')
#otree.Branch('b_flightSig3D', b_flightSig3D, 'b_flightSig3D/D')
#
#otree.Branch('jpsi_pt', jpsi_pt, 'jpsi_pt/D')
#otree.Branch('jpsi_eta', jpsi_eta, 'jpsi_eta/D')
#otree.Branch('jpsi_phi', jpsi_phi, 'jpsi_phi/D')
#otree.Branch('jpsi_mass', jpsi_mass, 'jpsi_mass/D')
#otree.Branch('jpsi_vprob', jpsi_vprob, 'jpsi_vprob/D')
#otree.Branch('jpsi_flightSig3D', jpsi_flightSig3D, 'jpsi_flightSig3D/D')

otree.Branch('tau_isRight_bS', tau_isRight_bS, 'tau_isRight_bS/O')
otree.Branch('tau_isRight_aS', tau_isRight_aS, 'tau_isRight_aS/O')
otree.Branch('tau_isRight_bS_ith', tau_isRight_bS_ith, 'tau_isRight_bS_ith/I')
otree.Branch('tau_isRight_aS_ith', tau_isRight_aS_ith, 'tau_isRight_aS_ith/I')
otree.Branch('tau_isRight_bS_n', tau_isRight_bS_n, 'tau_isRight_bS_n/I')
otree.Branch('tau_isRight_aS_n', tau_isRight_aS_n, 'tau_isRight_aS_n/I')

otree.Branch('event', event, 'event/I')
otree.Branch('run', run, 'run/I')
otree.Branch('lumi', lumi, 'lumi/I')

otree.Branch('vx', vx, 'vx/D')
otree.Branch('vy', vy, 'vy/D')
otree.Branch('vz', vz, 'vz/D')


Nevt = chain.GetEntries()

print 'Total Number of events = ', Nevt 
evtid = 0

for evt in xrange(Nevt):
    chain.GetEntry(evt)

    if evt%10000==0: print '{0:.2f}'.format(Double(evt)/Double(Nevt)*100.), '% processed'

#    import pdb; pdb.set_trace()


#    print len(tree.pft_tau_pt), tree.pft_tau_pt[0], tree.pft_tau_pt[1]

#    if len(chain.JpsiTau_ngentau)==0: continue
#    if chain.JpsiTau_ngentau[0]!=1: continue
    

#    print 'step2'

    if len(chain.JpsiTau_tau_pt)==0: 
#        print 'No reco. Tau ...'
        continue

#    print 'step3'
#    print 'saving ...'
    tau_pt[0] = chain.JpsiTau_tau_pt[0]
    tau_eta[0] = chain.JpsiTau_tau_eta[0]
    tau_phi[0] = chain.JpsiTau_tau_phi[0]
    tau_mass[0] = chain.JpsiTau_tau_mass[0]
    tau_vx[0] = chain.JpsiTau_tau_vx[0]
    tau_vy[0] = chain.JpsiTau_tau_vy[0]
    tau_vz[0] = chain.JpsiTau_tau_vz[0]
    tau_lip[0] = chain.JpsiTau_tau_lip[0]
    tau_lips[0] = chain.JpsiTau_tau_lips[0]
    tau_pvip[0] = chain.JpsiTau_tau_pvip[0]
    tau_pvips[0] = chain.JpsiTau_tau_pvips[0]
    tau_fl3d[0] = chain.JpsiTau_tau_fl3d[0]
    tau_fls3d[0] = chain.JpsiTau_tau_fls3d[0]
    tau_alpha[0] = chain.JpsiTau_tau_alpha[0]

    tau_max_dr_3prong[0] = chain.JpsiTau_tau_max_dr_3prong[0]
#        tau_iso03 [0] = chain.JpsiTau_tau_iso03[0]
    tau_charge[0] = chain.JpsiTau_tau_q[0]
    tau_vprob[0] = chain.JpsiTau_tau_vprob[0]
#        tau_dz2pv[0] = chain.JpsiTau_tau_dz2pv[0]
#        tau_dxy2pv[0] = chain.JpsiTau_tau_dxy2pv[0]
    tau_isRight[0] = chain.JpsiTau_tau_isRight[0]
    tau_isgen3matched[0] = int(bool(chain.JpsiTau_isgen3matched[0]))
    tau_isgen3[0] = int(bool(chain.JpsiTau_isgen3[0]))
#    tau_ngentau[0] = chain.JpsiTau_ngentau[0]
#        tau_nch[0] = chain.JpsiTau_nch[0]
#    else:
#        tau_pt[0] = -999
#        tau_eta[0] =  -999
#        tau_phi[0] =  -999
#        tau_mass[0] =  -999
#        tau_max_dr_3prong[0] =  -999
#        tau_iso03 [0] =  -999
#        tau_charge[0] =  -999
#        tau_flightSig3D[0] =  -999
#        tau_vprob[0] =  -999
#        tau_dz2pv[0] =  -999
#        tau_dxy2pv[0] =  -999
#        tau_isRight[0] =  -999
#        tau_isgen3matched[0] =  -999
#        tau_isgen3[0] =  -999
#        tau_ngentau[0] =  -999
#        tau_nch[0] =  -999
#
#
#    
    tau_gentaupt[0] = chain.JpsiTau_tau_matched_gentaupt[0]
#
#    if len(chain.JpsiTau_jpsi_pt)!=0 and len(chain.JpsiTau_tau_pt)!=0:
    b_mass[0] = chain.JpsiTau_B_mass[0]
    b_pt[0] =  chain.JpsiTau_B_pt[0]
    b_eta[0] =  chain.JpsiTau_B_eta[0]
    b_phi[0] =  chain.JpsiTau_B_phi[0]
#        b_flightSig3D[0] = chain.JpsiTau_b_flightSig3D[0]
#        
#        jpsi_pt[0] = chain.JpsiTau_jpsi_pt[0]
#        jpsi_eta[0] = chain.JpsiTau_jpsi_eta[0]
#        jpsi_phi[0] = chain.JpsiTau_jpsi_phi[0]
#        jpsi_mass[0] = chain.JpsiTau_jpsi_mass[0]
#        jpsi_vprob[0] = chain.JpsiTau_jpsi_vprob[0]
#        jpsi_flightSig3D[0] = chain.JpsiTau_jpsi_flightSig3D[0]
#    else:
#        
#        b_mass[0] = -999
#        b_pt[0] = -999
#        b_eta[0] = -999
#        b_phi[0] = -999
#        b_flightSig3D[0] = -999
#        
#        jpsi_pt[0] = -999
#        jpsi_eta[0] = -999
#        jpsi_phi[0] = -999
#        jpsi_mass[0] = -999
#        jpsi_vprob[0] = -999
#        jpsi_flightSig3D[0] = -999



#    if bool(chain.JpsiTau_tau_isRight[0]):
#        evtid += 1

#        print bool(chain.JpsiTau_tau_isRight[0])

    tau_isRight_bS[0] = chain.JpsiTau_isRight_bS[0]
    tau_isRight_aS[0] =  chain.JpsiTau_isRight_aS[0]
    tau_isRight_bS_ith[0] =  chain.JpsiTau_isRight_bS_ith[0]
    tau_isRight_aS_ith[0] =  chain.JpsiTau_isRight_aS_ith[0]
    tau_isRight_bS_n[0] =  chain.JpsiTau_isRight_bS_n[0]
    tau_isRight_aS_n[0] =  chain.JpsiTau_isRight_aS_n[0]

#    import pdb; pdb.set_trace()

    event[0] = chain.EVENT_event
    run[0] = chain.EVENT_run
    lumi[0] = chain.EVENT_lumiBlock

    vx[0] = chain.JpsiTau_bbPV_vx[0]
    vy[0] = chain.JpsiTau_bbPV_vy[0]
    vz[0] = chain.JpsiTau_bbPV_vz[0]

    evtid += 1
    otree.Fill()

otree.Write()
outputfile.Close()

print Nevt, 'evt processed.', evtid, 'evt has matching'

