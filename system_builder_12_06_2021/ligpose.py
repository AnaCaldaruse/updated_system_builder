#!/usr/bin/env python
# coding: utf-8
#
# The main change that is needed here is that it needs to:
#    read in the charged molecule from ligand_builder
# ----------------------------------------------------------------------
# ligpose.py, version 1.0, Mary Pitman, Mobley Lab UCI
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appears in supporting documentation.
# ----------------------------------------------------------------------
__doc__="""
        Modules designed to posit ligands using the OEPosit space of
        Openeye. This module performs ligand posing through
        information from existing complex structures or through
        docking.
        
        arguments for main(): (pdb_path, mtz_path, lig_in, method)
        Input arguments are str. Example: 'shapefit'
        
        pdb_path: path to the *.pdb file containing your reference
                  complex structure.
        mtz_path: *.mtz file for the reference complex structure.
                  *.mtz files can be downloaded from the PDB.
        lig_in:   The *.mol2 file generated from ligbuild.py or
                  through other means. Can take *.pdb and *.sdf.
                  Any 3D structure.
        method:   the method you wish to use for OEPosit().
                  shapefit, mcs, hybrid, fred."
                
        
        Outputs:
            - *.sdf and *.mol2 files for building simulations.
            - *.svg files to visualize posed ligand interactions and
              compare to the reference complex interactions.
            - a receptor.pdb to generate the updated bound complex.
        """
__version__ = '1.0'
__author__ = 'M. Pitman'

import os
import sys
import re
import subprocess as sub

from openeye import oechem
from openeye import oespruce
from openeye import oedocking
from openeye import oegrid

from oescripts import du2liginters
toolkit_root = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------

def extract_ligand(du, ofs):
    ''' Extracts and checks for a ligand in a complex.'''
    ligand = oechem.OEGraphMol()
    if not du.HasLigand():
        oechem.OEThrow.Fatal("Error: There is no ligand in the"
                             " OEDesignUnit."
                             )
    oechem.OEWriteMolecule(ofs, ligand)
    ofs.close()
    
    return ligand
    

def build_du(pdb_path):
#def build_du(pdb_path, mtz_path):
    ''' Builds the initial design unit. '''
    # Read the input pdb file 
    iname = os.path.abspath(pdb_path)
    
    #read PDB via oemolistream (supports SetFlavor)
    ims = oechem.oemolistream()
    if not ims.open(iname):
        oechem.OEThrow.Fatal(f"Unable to open {iname} for reading")

    ims.SetFlavor(
        oechem.OEFormat_PDB,
        oechem.OEIFlavor_PDB_Default
        | oechem.OEIFlavor_PDB_DATA
        | oechem.OEIFlavor_PDB_ALTLOC
    )

    complexmol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ims, complexmol):
        oechem.OEThrow.Fatal(f"Unable to read {iname}")

    ims.close()

    # Spruce DU construction
    metadata = oespruce.OEStructureMetadata()
    opts = oespruce.OEMakeDesignUnitOptions()
    dus = oespruce.OEMakeDesignUnits(complexmol, metadata, opts)
    
    # pick first DU
    du = None
    for du_candidate in dus:
        du = du_candidate
        break

    if du is None:
        oechem.OEThrow.Fatal("No design units were created")
    
    # protonate DU
    oespruce.OEProtonateDesignUnit(du)

    du_ofile = os.path.basename(iname)[:-4] + "_DU.oedu"
    ofs = oechem.oeofstream(du_ofile)
    oechem.OEWriteDesignUnit(ofs, du)
    ofs.close()

    print(f"\n--------------------------------------------------\n"
          f"Saved a design unit binary file to {du_ofile}")
    
    return du, du_ofile
 

def make_receptor(oedu_path):
    '''
    Defines how to make the receptor that will be
    added to the design unit.
    '''
    recOpts = oedocking.OEMakeReceptorOptions()
    iname = oedu_path
    ifs = oechem.oeifstream()
    ifs.open(iname)

    rec_ofile = 'receptor_out.oedu'
    ofs = oechem.oeofstream(rec_ofile)

    du = oechem.OEDesignUnit()
    while oechem.OEReadDesignUnit(ifs, du):
        if oedocking.OEMakeReceptor(du, recOpts):
            oechem.OEWriteDesignUnit(ofs, du)
        else:
            oechem.OEThrow.Warning("Failed to make receptor.")
    
    ifs.close()
    ofs.close()
    
    return rec_ofile
    

def flexible_overlay(pdb_path, lig_in, method):
    '''
    Positions the ligand, sets options for OEPosit(),
    and outputs ligand and receptor structure files.
    Returns:
        posed .sdf path on success
        None on failure (instead of crashing)
    '''

    # Build design unit from receptor PDB
    du, du_ofile = build_du(pdb_path)
    rec_ofile = make_receptor(du_ofile)

    # Read output receptor
    ifs = oechem.oeifstream()
    if not ifs.open(rec_ofile):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % rec_ofile)

    du = oechem.OEDesignUnit()
    oechem.OEReadDesignUnit(ifs, du)
    ifs.close()

    if du.HasReceptor():
        print("Receptor built for the design unit, du.")
    else:
        print("The design unit does not have a receptor built")

    # Read in the ligand to be positioned
    try:
        lig_to_pose = oechem.OEMol(lig_in)
    except TypeError:
        ims = oechem.oemolistream(lig_in)
        lig_to_pose = oechem.OEMol()
        oechem.OEReadMolecule(ims, lig_to_pose)
        ims.close()

    if lig_to_pose.GetTitle() == '':
        lig_to_pose.SetTitle('posed_lig')
        print("The ligand to be posed did not have a title. "
              "Ligand title set to posed_lig.")

    print("Posing in progress ...")

    # First, try shapefit
    opts = oedocking.OEPositOptions()
    opts.SetFullConformationSearch(True)
    opts.SetIgnoreNitrogenStereo(True)
    opts.SetPositMethods(oedocking.OEPositMethod_SHAPEFIT)

    poser = oedocking.OEPosit(opts)
    poser.AddReceptor(du)

    result = oedocking.OESinglePoseResult()
    code = poser.Dock(result, lig_to_pose)

    if code == oedocking.OEDockingReturnCode_Success:
        print("[pose] shapefit success")
    else:
        print("[pose] shapefit failed: "
              f"{oedocking.OEDockingReturnCodeGetName(code)}")

        # If shapefit fails, try hybrid 
        opts = oedocking.OEPositOptions()
        opts.SetFullConformationSearch(True)
        opts.SetIgnoreNitrogenStereo(True)
        opts.SetPositMethods(oedocking.OEPositMethod_HYBRID)

        poser = oedocking.OEPosit(opts)
        poser.AddReceptor(du)

        result = oedocking.OESinglePoseResult()
        code = poser.Dock(result, lig_to_pose)

        if code == oedocking.OEDockingReturnCode_Success:
            print("[pose] hybrid fallback success")
        else:
            print("[pose] hybrid fallback failed: "
                  f"{oedocking.OEDockingReturnCodeGetName(code)}")

            # Lastly, if hybrid also fails try fred for more traditional docking
            opts = oedocking.OEPositOptions()
            opts.SetFullConformationSearch(True)
            opts.SetIgnoreNitrogenStereo(True)
            opts.SetPositMethods(oedocking.OEPositMethod_FRED)

            poser = oedocking.OEPosit(opts)
            poser.AddReceptor(du)

            result = oedocking.OESinglePoseResult()
            code = poser.Dock(result, lig_to_pose)

            if code == oedocking.OEDockingReturnCode_Success:
                print("[pose] conformer fallback success")
            else:
                print("[pose] conformer fallback failed: "
                      f"{oedocking.OEDockingReturnCodeGetName(code)}")
                print(f"[pose] total failure → no poses for {lig_to_pose.GetTitle()}")
                return None

    # Success?
    posed_du = result.GetDesignUnit()
    posed_oemol = oechem.OEGraphMol()
    posed_du.GetLigand(posed_oemol)
    du_protein = oechem.OEGraphMol()
    posed_du.GetProtein(du_protein)

    title = lig_to_pose.GetTitle()

    # Save posed ligand
    for ext in ("sdf", "mol2"):
        ofile = f"posed_{title}.{ext}"
        ostream = oechem.oemolostream()
        ostream.open(ofile)
        oechem.OEWriteMolecule(ostream, posed_oemol)
        ostream.close()

    # Save receptor (optional for visualization/debug)
    ostream = oechem.oemolostream("receptor_protein.pdb")
    oechem.OEWriteMolecule(ostream, du_protein)
    ostream.close()

    print(f"[pose] final success for {title}")

    # Optional: interaction fingerprints
    try:
        du2liginters(posed_du)
    except Exception:
        pass

    return f"posed_{title}.sdf"

  
def visualize_contacts(pdb_path, lig_out):
    '''
    Uses an openeye script contained in oescripts to create
    a visualization of how the ligand interacts with the structure
    in the reference complex and the generated posed complex.
    
    example of how to open:
        open -a "Google Chrome" ouput.svg
    '''
    # Create the visualization of the reference complex
    # interactions.
    resolved_inter_ofile = os.path.basename(pdb_path)[:-4] \
                           + "_interactions.svg"
    script = os.path.join(toolkit_root, "oescripts", "complex2img.py")
    sub.call(f"python3 {script} -complex {pdb_path} -out {resolved_inter_ofile}", shell=True)

    print("\n--------------------------------------------------\n"
          "The interactions between the reference complex "
          "structure and resolved ligand were visualized at:\n"
          f"{resolved_inter_ofile}\n"
    )

    # Create the visualization of the posed ligand interactions.
    ofile2 = 'receptor_protein.pdb'
    base_ofile = lig_out.split(".")[0]
    inter_ofile = f"{base_ofile}_interactions.svg"
    sub.call(
        f"python3 {script} -protein {ofile2} -ligand {lig_out} -out {inter_ofile}",
        shell=True
        )
    print("The interactions between the posed ligand and "
          "the receptor were visualized at:\n"
          f"{inter_ofile}\n"
          "--------------------------------------------------"
    )


def clean_up(files):
    '''
    Removes *.oedu files generated by ligpose.py.
    Input: list of files to remove.
    '''
    print("Cleaning up files from ligand posing...")
    for i in files:
        if os.path.exists(i):
            os.remove(i)
            print("Design unit binary file {} removed by"
                  " ligpose.clean_up".format(i))
        else:
            print("The file {} does not exist for clean up."
                  .format(i))

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("ligpose.log", "a", buffering=1)  # line-buffered

    def write(self, outputs):
        self.terminal.write(outputs)
        self.log.write(outputs)
        self.log.flush()   # <-- critical!

    def flush(self):
        self.log.flush()

 
#def main(pdb_path, mtz_path, lig_in, method, generate_svgs=True):
def main(pdb_path, lig_in, method, generate_svgs=True):
    lig_out = flexible_overlay(pdb_path, lig_in, method)
    if lig_out is None:
        return None

    if generate_svgs:
        visualize_contacts(pdb_path, lig_out)
    # Remove the generated binary files.
    to_clean = [os.path.basename(pdb_path)[:-4] + "_DU.oedu",
                'receptor_out.oedu'
                ]
    clean_up(to_clean)
    
    return lig_out

if __name__ == '__main__':
    main()
