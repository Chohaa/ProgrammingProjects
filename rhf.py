  1 from functools import reduce
Last login: Mon Mar 20 13:39:13 on console
xcrun: error: invalid active developer path (/Library/Developer/CommandLineTools), missing xcrun at: /Library/Developer/CommandLineTools/usr/bin/xcrun
 hanna@macarena2  ~
 vi fermi_rhf.py
 ✘ hanna@macarena2  ~
 cat fermi_rhf.py
from functools import reduce
import numpy as np
import scipy

import pyscf
from pyscf import fci
from pyscf import gto, scf, ao2mo, lo, tdscf, cc


def tda_denisty_matrix(td, state_id):
    '''
    Taking the TDA amplitudes as the CIS coefficients, calculate the density
    matrix (in AO basis) of the excited states
    '''
    cis_t1 = td.xy[state_id][0]
    dm_oo =-np.einsum('ia,ka->ik', cis_t1.conj(), cis_t1)
    dm_vv = np.einsum('ia,ic->ac', cis_t1, cis_t1.conj())

    # The ground state density matrix in mo_basis
    mf = td._scf
    dm = np.diag(mf.mo_occ)

    # Add CIS contribution
    nocc = cis_t1.shape[0]
    # Note that dm_oo and dm_vv correspond to spin-up contribution. "*2" to
    # include the spin-down contribution
    dm[:nocc,:nocc] += dm_oo * 2
    dm[nocc:,nocc:] += dm_vv * 2

    # Transform density matrix to AO basis
    mo = mf.mo_coeff
    dm = np.einsum('pi,ij,qj->pq', mo, dm, mo.conj())
    return dm

mol = gto.Mole()
mol.atom = '''
H       -3.4261000000     -2.2404000000      5.4884000000
H       -5.6274000000     -1.0770000000      5.2147000000
C       -3.6535000000     -1.7327000000      4.5516000000
H       -1.7671000000     -2.2370000000      3.6639000000
C       -4.9073000000     -1.0688000000      4.3947000000
H       -6.1631000000      0.0964000000      3.1014000000
C       -2.7258000000     -1.7321000000      3.5406000000
H       -0.3003000000      1.0832000000     -5.2357000000
C       -5.2098000000     -0.4190000000      3.2249000000
C       -2.9961000000     -1.0636000000      2.3073000000
H       -1.1030000000     -1.5329000000      1.3977000000
H       -0.4270000000     -0.8029000000     -0.8566000000
H        0.2361000000     -0.0979000000     -3.1273000000
C       -1.0193000000      1.0730000000     -4.4150000000
H       -2.4988000000      2.2519000000     -5.5034000000
C       -4.2740000000     -0.3924000000      2.1445000000
H       -5.5015000000      0.7944000000      0.8310000000
C       -2.0613000000     -1.0272000000      1.2718000000
C       -1.3820000000     -0.2895000000     -0.9772000000
C       -0.7171000000      0.4180000000     -3.2476000000
C       -2.2720000000      1.7395000000     -4.5690000000
H       -4.1576000000      2.2412000000     -3.6787000000
C       -4.5463000000      0.2817000000      0.9534000000
C       -2.3243000000     -0.3402000000      0.0704000000
 59 n_singlets = 1
 70     # The ground state density matrix in mo_basis
C       -1.6528000000      0.3874000000     -2.1670000000
C       -3.1998000000      1.7341000000     -3.5584000000
C       -3.6044000000      0.3309000000     -0.0943000000
C       -2.9302000000      1.0591000000     -2.3292000000
C       -3.8665000000      1.0187000000     -1.2955000000
H       -4.8243000000      1.5256000000     -1.4217000000
'''

mol.basis = '6-31g*'
mol.spin = 0
mol.build()

mf = scf.RHF(mol)
mf.verbose = 4
mf.get_init_guess(mol, key='minao')
mf.conv_tol = 1e-9
mf.level_shift = .1
mf.diis_start_cycle = 4
mf.diis_space = 10
mf.run(max_cycle=200)
#mf.chkfile.dump_scf(mol, 'tet-4mer', mf.e_tot, mf.mo_energy, mf.mo_coeff, mf.mo_occ, overwrite_mol=True)

n_triplets = 2
n_singlets = 1

avg_rdm1 = mf.make_rdm1()


# compute singlets
mytd = tdscf.TDA(mf)
mytd.singlet = True
mytd = mytd.run(nstates=n_singlets)
mytd.analyze()
for i in range(mytd.nroots):
    avg_rdm1 += tda_denisty_matrix(mytd, i)

# compute triplets
mytd = tdscf.TDA(mf)
mytd.singlet = False
mytd = mytd.run(nstates=n_triplets)
mytd.analyze()
for i in range(mytd.nroots):
    avg_rdm1 += tda_denisty_matrix(mytd, i)

# normalize
avg_rdm1 = avg_rdm1 / (n_singlets + n_triplets + 1)


S = mf.get_ovlp()
np.save("overlap_mat", S)
np.save("density_mat", mf.make_rdm1())
np.save("rhf_mo_coeffs", mf.mo_coeff)
np.save("cis_sa_density_mat", avg_rdm1)


 hanna@macarena2  ~
 vi fermi_rhf.py
 hanna@macarena2  ~
 vi fermi_rhf.py
 hanna@macarena2  ~
 ls
Class                  Downloads              Music                  SingletFission         iCloud Drive (Archive) pyscf
Computer               Julia                  Pictures               TetraceneActiveSpace   oh-my-zsh
Desktop                Library                ProgrammingProjects    Zotero                 opt
Documents              Movies                 Public                 get-pip.py             pauli_jupyter.ipynb
 hanna@macarena2  ~
 cd TetraceneActiveSpace
 hanna@macarena2  ~/TetraceneActiveSpace
 cd HF
 hanna@macarena2  ~/TetraceneActiveSpace/HF
 ls
cis_sa_density_mat.npy density_mat.npy        fermi_rhf.out          fermi_rhf.py           overlap_mat.npy        rhf_mo_coeffs.npy
 hanna@macarena2  ~/TetraceneActiveSpace/HF
 vi fermi_rhf.py
 ✘ hanna@macarena2  ~/TetraceneActiveSpace/HF
 vi fermi_rhf.py
 hanna@macarena2  ~/TetraceneActiveSpace/HF
 cat fermi_rhf.py
from functools import reduce
import numpy as np
import scipy

import pyscf
from pyscf import fci
from pyscf import gto, scf, ao2mo, lo, tdscf, cc


def tda_denisty_matrix(td, state_id):
    '''
    Taking the TDA amplitudes as the CIS coefficients, calculate the density
    matrix (in AO basis) of the excited states
    '''
    cis_t1 = td.xy[state_id][0]
    dm_oo =-np.einsum('ia,ka->ik', cis_t1.conj(), cis_t1)
    dm_vv = np.einsum('ia,ic->ac', cis_t1, cis_t1.conj())

    # The ground state density matrix in mo_basis
    mf = td._scf
    dm = np.diag(mf.mo_occ)

    # Add CIS contribution
    nocc = cis_t1.shape[0]
    # Note that dm_oo and dm_vv correspond to spin-up contribution. "*2" to
    # include the spin-down contribution
    dm[:nocc,:nocc] += dm_oo * 2
    dm[nocc:,nocc:] += dm_vv * 2

    # Transform density matrix to AO basis
    mo = mf.mo_coeff
    dm = np.einsum('pi,ij,qj->pq', mo, dm, mo.conj())
    return dm

mol = gto.Mole()
mol.atom = '''
H       -3.4261000000     -2.2404000000      5.4884000000
H       -5.6274000000     -1.0770000000      5.2147000000
C       -3.6535000000     -1.7327000000      4.5516000000
H       -1.7671000000     -2.2370000000      3.6639000000
C       -4.9073000000     -1.0688000000      4.3947000000
H       -6.1631000000      0.0964000000      3.1014000000
C       -2.7258000000     -1.7321000000      3.5406000000
H       -0.3003000000      1.0832000000     -5.2357000000
C       -5.2098000000     -0.4190000000      3.2249000000
C       -2.9961000000     -1.0636000000      2.3073000000
H       -1.1030000000     -1.5329000000      1.3977000000
H       -0.4270000000     -0.8029000000     -0.8566000000
H        0.2361000000     -0.0979000000     -3.1273000000
C       -1.0193000000      1.0730000000     -4.4150000000
H       -2.4988000000      2.2519000000     -5.5034000000
C       -4.2740000000     -0.3924000000      2.1445000000
H       -5.5015000000      0.7944000000      0.8310000000
C       -2.0613000000     -1.0272000000      1.2718000000
C       -1.3820000000     -0.2895000000     -0.9772000000
C       -0.7171000000      0.4180000000     -3.2476000000
C       -2.2720000000      1.7395000000     -4.5690000000
H       -4.1576000000      2.2412000000     -3.6787000000
C       -4.5463000000      0.2817000000      0.9534000000
C       -2.3243000000     -0.3402000000      0.0704000000
C       -1.6528000000      0.3874000000     -2.1670000000
C       -3.1998000000      1.7341000000     -3.5584000000
C       -3.6044000000      0.3309000000     -0.0943000000
C       -2.9302000000      1.0591000000     -2.3292000000
C       -3.8665000000      1.0187000000     -1.2955000000
H       -4.8243000000      1.5256000000     -1.4217000000
'''

mol.basis = '6-31g*'
mol.spin = 0
mol.build()

mf = scf.RHF(mol)
mf.verbose = 4
mf.get_init_guess(mol, key='minao')
mf.conv_tol = 1e-9
mf.level_shift = .1
mf.diis_start_cycle = 4
mf.diis_space = 10
mf.run(max_cycle=200)
#mf.chkfile.dump_scf(mol, 'tet-4mer', mf.e_tot, mf.mo_energy, mf.mo_coeff, mf.mo_occ, overwrite_mol=True)

n_triplets = 2
n_singlets = 1

avg_rdm1 = mf.make_rdm1()


# compute singlets
mytd = tdscf.TDA(mf)
mytd.singlet = True
mytd = mytd.run(nstates=n_singlets)
mytd.analyze()
for i in range(mytd.nroots):
    avg_rdm1 += tda_denisty_matrix(mytd, i)

  1 from pyscf import gto, scf, dft, tddft
# compute triplets
mytd = tdscf.TDA(mf)
mytd.singlet = False
mytd = mytd.run(nstates=n_triplets)
mytd.analyze()
for i in range(mytd.nroots):
    avg_rdm1 += tda_denisty_matrix(mytd, i)
  1 from pyscf import gto, scf, dft, tddft

# normalize
avg_rdm1 = avg_rdm1 / (n_singlets + n_triplets + 1)


S = mf.get_ovlp()
np.save("overlap_mat", S)
np.save("density_mat", mf.make_rdm1())
np.save("rhf_mo_coeffs", mf.mo_coeff)
np.save("cis_sa_density_mat", avg_rdm1)


 hanna@macarena2  ~/TetraceneActiveSpace/HF
 ls
2,2.out                cis_sa_density_mat.npy fermi_rhf.py           integrals_h2.npy       overlap_mat.npy
2,2.py                 density_mat.npy        integrals_h0.npy       mo_coeffs_act.npy      rhf_mo_coeffs.npy
Cact.molden            fermi_rhf.out          integrals_h1.npy       mo_coeffs_doc.npy      xyz.npy
 hanna@macarena2  ~/TetraceneActiveSpace/HF
 python
mytd = tddft.TDDFT(mf, frozen=[0:58]).run() >tddft.out
Python 3.9.13 (main, Aug 25 2022, 18:29:29)
[Clang 12.0.0 ] :: Anaconda, Inc. on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> exit()
  1 from functools import reduce
zsh: no matches found: tddft.TDDFT(mf, frozen=[0:58]).run
 ✘ hanna@macarena2  ~/TetraceneActiveSpace/HF
 python mytd = tddft.TDDFT(mf, frozen=[0:58]).run() >tddft.out
zsh: no matches found: tddft.TDDFT(mf, frozen=[0:58]).run
 ✘ hanna@macarena2  ~/TetraceneActiveSpace/HF
 python
Python 3.9.13 (main, Aug 25 2022, 18:29:29)
[Clang 12.0.0 ] :: Anaconda, Inc. on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> mytd = tddft.TDDFT(mf, frozen=[0:58]).run() >tddft.out
  File "<stdin>", line 1
    mytd = tddft.TDDFT(mf, frozen=[0:58]).run() >tddft.out
                                    ^
SyntaxError: invalid syntax
>>> ls
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'ls' is not defined
>>> pwd
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'pwd' is not defined
>>> ;
  File "<stdin>", line 1
    ;
    ^
SyntaxError: invalid syntax
>>> [
... pwd
... ls
  File "<stdin>", line 3
    ls
    ^
SyntaxError: invalid syntax
>>> exit()
 hanna@macarena2  ~/TetraceneActiveSpace/HF
 from pyscf import gto, scf, dft, tddft
zsh: command not found: from
 ✘ hanna@macarena2  ~/TetraceneActiveSpace/HF
 python from pyscf import gto, scf, dft, tddft
python: can't open file '/Users/hanna/TetraceneActiveSpace/HF/from': [Errno 2] No such file or directory
 ✘ hanna@macarena2  ~/TetraceneActiveSpace/HF
 vi dft
 hanna@macarena2  ~/TetraceneActiveSpace/HF
 ls
2,2.out                cis_sa_density_mat.npy fermi_rhf.out          integrals_h1.npy       mo_coeffs_doc.npy      xyz.npy
2,2.py                 density_mat.npy        fermi_rhf.py           integrals_h2.npy       overlap_mat.npy
Cact.molden            dft                    integrals_h0.npy       mo_coeffs_act.npy      rhf_mo_coeffs.npy
 hanna@macarena2  ~/TetraceneActiveSpace/HF
 vi dft
 hanna@macarena2  ~/TetraceneActiveSpace/HF
 pythong dft.py
zsh: command not found: pythong
 ✘ hanna@macarena2  ~/TetraceneActiveSpace/HF
 python dft.py
  File "/Users/hanna/TetraceneActiveSpace/HF/dft.py", line 3
    mytd = tddft.TDDFT(mf, frozen=[0:58]).run()
                                    ^
SyntaxError: invalid syntax
 ✘ hanna@macarena2  ~/TetraceneActiveSpace/HF
 ls
2,2.out                cis_sa_density_mat.npy dft.py                 integrals_h0.npy       mo_coeffs_act.npy      rhf_mo_coeffs.npy
2,2.py                 density_mat.npy        fermi_rhf.out          integrals_h1.npy       mo_coeffs_doc.npy      xyz.npy
Cact.molden            dft                    fermi_rhf.py           integrals_h2.npy       overlap_mat.npy
 hanna@macarena2  ~/TetraceneActiveSpace/HF
 rm dft
 hanna@macarena2  ~/TetraceneActiveSpace/HF
 python dft.py
Traceback (most recent call last):
  File "/Users/hanna/TetraceneActiveSpace/HF/dft.py", line 3, in <module>
    mytd = tddft.TDDFT(mf, frozen=[58]).run()
NameError: name 'mf' is not defined
 ✘ hanna@macarena2  ~/TetraceneActiveSpace/HF
 vi 2,2.py
 hanna@macarena2  ~/TetraceneActiveSpace/HF
 cat 2,2.py
from functools import reduce
import numpy as np
import scipy

import pyscf
from pyscf import fci
from pyscf import gto, scf, ao2mo, lo, tdscf, cc, dft, tddft



def get_natural_orbital_active_space(rdm, S, thresh=.01):


    Ssqrt = scipy.linalg.sqrtm((S+S.T)/2.0)
    Sinvsqrt = scipy.linalg.inv(Ssqrt)

    print(" Number of electrons found %12.8f" %np.trace(S@rdm))

    Dtot = Ssqrt.T @ rdm @ Ssqrt
    #Dtot = Ssqrt.T @ ( da + db) @ Ssqrt
    D_evals, D_evecs = np.linalg.eigh((Dtot+Dtot.T)/2.0)

    sorted_list = np.argsort(D_evals)[::-1]
    D_evals = D_evals[sorted_list]
    D_evecs = D_evecs[:,sorted_list]

    act_list = []
    doc_list = []


    for idx,n in enumerate(D_evals):
        print(" %4i = %12.8f" %(idx,n),end="")
        if n < 2.0 - thresh:
            if n > thresh:
                act_list.append(idx)
                print(" Active")
            else:
                print(" Virt")
        else:
            doc_list.append(idx)
            print(" DOcc")

    print(" Number of active orbitals: ", len(act_list))
    print(" Number of doc    orbitals: ", len(doc_list))

    D_evecs = Sinvsqrt @ D_evecs
    Cdoc = D_evecs[:, doc_list]
    Cact = D_evecs[:, act_list]
    return Cdoc, Cact


mol = gto.Mole()
mol.atom = '''
H       -3.4261000000     -2.2404000000      5.4884000000
H       -5.6274000000     -1.0770000000      5.2147000000
C       -3.6535000000     -1.7327000000      4.5516000000
H       -1.7671000000     -2.2370000000      3.6639000000
C       -4.9073000000     -1.0688000000      4.3947000000
H       -6.1631000000      0.0964000000      3.1014000000
C       -2.7258000000     -1.7321000000      3.5406000000
H       -0.3003000000      1.0832000000     -5.2357000000
C       -5.2098000000     -0.4190000000      3.2249000000
C       -2.9961000000     -1.0636000000      2.3073000000
H       -1.1030000000     -1.5329000000      1.3977000000
H       -0.4270000000     -0.8029000000     -0.8566000000
H        0.2361000000     -0.0979000000     -3.1273000000
C       -1.0193000000      1.0730000000     -4.4150000000
H       -2.4988000000      2.2519000000     -5.5034000000
C       -4.2740000000     -0.3924000000      2.1445000000
H       -5.5015000000      0.7944000000      0.8310000000
C       -2.0613000000     -1.0272000000      1.2718000000
C       -1.3820000000     -0.2895000000     -0.9772000000
C       -0.7171000000      0.4180000000     -3.2476000000
C       -2.2720000000      1.7395000000     -4.5690000000
H       -4.1576000000      2.2412000000     -3.6787000000
C       -4.5463000000      0.2817000000      0.9534000000
C       -2.3243000000     -0.3402000000      0.0704000000
126 n_triplets = 2
C       -1.6528000000      0.3874000000     -2.1670000000
C       -3.1998000000      1.7341000000     -3.5584000000
 97     D_evals = D_evals[sorted_list]
C       -3.6044000000      0.3309000000     -0.0943000000
C       -2.9302000000      1.0591000000     -2.3292000000
109                 print(" Active")
C       -3.8665000000      1.0187000000     -1.2955000000
H       -4.8243000000      1.5256000000     -1.4217000000
113             doc_list.append(idx)
105         print(" %4i = %12.8f" %(idx,n),end="")
'''

  1 from functools import reduce

np.save("xyz.npy", mol.atom)

mol.basis = '6-31g*'
mol.spin = 0
mol.build()

mf = scf.RHF(mol).density_fit()
  1 from functools import reduce
mf.verbose = 4
mf.get_init_guess(mol, key='minao')
mf.conv_tol = 1e-9

115     print(" Number of doc    orbitals: ", len(doc_list))
# load precomputed data
C = np.load("rhf_mo_coeffs.npy")
121
avg_rdm1 = np.load("cis_sa_density_mat.npy")

170 # Rotate 1electron terms to active space
171 h = Cact.T @ h @ Cact
172 j = Cact.T @ j @ Cact;
173 k = Cact.T @ k @ Cact;
174
175 h1 = h + j - .5*k;
176
177 # form 2e integrals in active space
178 nact = h.shape[0]
179 h2 = pyscf.ao2mo.kernel(mol, Cact, aosym="s4", compact=False)
180 h2.shape = (nact, nact, nact, nact)
181
182 #mytd = tddft.TDDFT(mf, frozen=2).run()
183
184 np.save("integrals_h0", h0)
185 np.save("integrals_h1", h1)
186 np.save("integrals_h2", h2)
187 np.save("mo_coeffs_act", Cact)
188 np.save("mo_coeffs_doc", Cdoc)
189
190 np.save("overlap_mat", S)
191 np.save("density_mat", avg_rdm1)
192 np.save("rhf_mo_coeffs", mf.mo_coeff)
193
194
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~                                                                                             ~/TetraceneActiveSpace/HF/fermi_rhf.py\
