

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
