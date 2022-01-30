import os
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign


def cs_confgen(input, output, numconf, add_ref):
    param = AllChem.ETKDGv3()
    param.pruneRmsThresh = 0.5
    param.useSmallRingTorsions = True
    param.numThreads=0

    name, ext = os.path.splitext(os.path.basename(input))
    if ext.lower() in ['.smi']:
        mols = Chem.SmilesMolSupplier(input, " \t", 0, 1, titleLine=True, sanitize=True)
        print("SMILES are being processing")
    elif ext.lower() in ['.sdf', 'sd']:
        mols = Chem.SDMolSupplier(input, removeHs=True)
        print("SDF are being processing")
    else:
        print("Uknown extention")
        exit()

    for mol in mols:
        nr = int(AllChem.CalcNumRotatableBonds(mol))
        if nr <= 3:
            nc = 50
        elif nr > 6:
            nc = 300
        else:
            nc = nr**3

        if numconf:
            nc = numconf
        label = mol.GetProp("_Name")
        print("NC", nc, "of ", label)

        mol = Chem.AddHs(mol)
        cids = AllChem.EmbedMultipleConfs(mol, numConfs = nc, params = param)
        mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=8192, mmffVariant='MMFF94')

        sdfile = open(output, 'at')
        w = Chem.SDWriter(sdfile)
        res = []
        for cid in cids:
            ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
            e = ff.CalcEnergy()
            res.append((cid, e))
        sorted_res = sorted(res, key=lambda x:x[1])
        rdMolAlign.AlignMolConformers(mol)
        for cid, e in sorted_res:
            mol.SetProp('_Name', str(label))
            mol.SetProp('idnumber', str(label))
            mol.SetProp('CID', str(cid))
            mol.SetProp('Energy', str(e))
            w.write(mol, confId=cid)
        w.close()

if __name__=='__main__':
    cs_confgen()
