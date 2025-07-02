# python LHEReader_v1p0.py --inputs file1.lhe file2.lhe file3.lhe --output combined.root

# LHE READER - Multiple File Support
import xml.etree.ElementTree as ET
import ROOT
from ROOT import TFile, TTree, std
from array import array
import argparse
import pathlib

def read_xml(filename):
    tree = ET.parse(filename)
    return tree.getroot()

def read_xml_child(root):
    return [child.text.split() for child in root if child.tag == "event"]

def setup_tree():
    tree = TTree("events", "events")
    m_Npart = array("i", [0])
    m_eventweight = array("f", [0.0])
    m_scale = array("f", [0.0])
    m_qed = array("f", [0.0])
    m_qcd = array("f", [0.0])
    m_pid = std.vector("int")()
    m_status = std.vector("int")()
    m_mother1 = std.vector("int")()
    m_mother2 = std.vector("int")()
    m_color1 = std.vector("int")()
    m_color2 = std.vector("int")()
    m_px = std.vector("float")()
    m_py = std.vector("float")()
    m_pz = std.vector("float")()
    m_e = std.vector("float")()
    m_mass = std.vector("float")()
    m_tau = std.vector("float")()
    m_spin = std.vector("float")()

    tree.Branch("numParticles", m_Npart, "numParticles/I")
    tree.Branch("eventweight", m_eventweight, "eventweight/F")
    tree.Branch("scale", m_scale, "scale/F")
    tree.Branch("alpha_qed", m_qed, "alpha_qed/F")
    tree.Branch("alpha_qcd", m_qcd, "alpha_qcd/F")
    tree.Branch("pid", m_pid)
    tree.Branch("status", m_status)
    tree.Branch("mother1", m_mother1)
    tree.Branch("mother2", m_mother2)
    tree.Branch("color1", m_color1)
    tree.Branch("color2", m_color2)
    tree.Branch("px", m_px)
    tree.Branch("py", m_py)
    tree.Branch("pz", m_pz)
    tree.Branch("energy", m_e)
    tree.Branch("mass", m_mass)
    tree.Branch("tau", m_tau)
    tree.Branch("spin", m_spin)

    branches = (m_Npart, m_eventweight, m_scale, m_qed, m_qcd,
                m_pid, m_status, m_mother1, m_mother2, m_color1,
                m_color2, m_px, m_py, m_pz, m_e, m_mass, m_tau, m_spin)
    return tree, branches

def fill_tree(data, tree, branches):
    (m_Npart, m_eventweight, m_scale, m_qed, m_qcd,
     m_pid, m_status, m_mother1, m_mother2, m_color1,
     m_color2, m_px, m_py, m_pz, m_e, m_mass, m_tau, m_spin) = branches

    for event in data:
        m_Npart[0] = int(event[0])
        m_eventweight[0] = float(event[2])
        m_scale[0] = float(event[3])
        m_qed[0] = float(event[4])
        m_qcd[0] = float(event[5])
        x, y = 6, 19

        for _ in range(m_Npart[0]):
            m_pid.push_back(int(event[x]))
            m_status.push_back(int(event[x + 1]))
            m_mother1.push_back(int(event[x + 2]))
            m_mother2.push_back(int(event[x + 3]))
            m_color1.push_back(int(event[x + 4]))
            m_color2.push_back(int(event[x + 5]))
            m_px.push_back(float(event[x + 6]))
            m_py.push_back(float(event[x + 7]))
            m_pz.push_back(float(event[x + 8]))
            m_e.push_back(float(event[x + 9]))
            m_mass.push_back(float(event[x + 10]))
            m_tau.push_back(float(event[x + 11]))
            m_spin.push_back(float(event[x + 12]))
            x, y = y, y + 13

        tree.Fill()
        m_pid.clear(); m_status.clear(); m_mother1.clear(); m_mother2.clear()
        m_color1.clear(); m_color2.clear(); m_px.clear(); m_py.clear(); m_pz.clear()
        m_e.clear(); m_mass.clear(); m_tau.clear(); m_spin.clear()

def main(input_files, output_file):
    tree, branches = setup_tree()

    for input_file in input_files:
        print(f"Processing: {input_file}")
        root = read_xml(input_file)
        data = read_xml_child(root)
        fill_tree(data, tree, branches)

    with TFile.Open(output_file, "recreate") as f:
        tree.Write("", TFile.kOverwrite)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge multiple LHE files into a ROOT TTree")
    parser.add_argument("--inputs", nargs='+', required=True, help="List of input .lhe files")
    parser.add_argument("--output", required=True, help="Output .root file")

    args = parser.parse_args()
    input_files = [str(path) for path in args.inputs if pathlib.Path(path).suffix == ".lhe"]

    if not input_files:
        print("No valid .lhe files provided.")
    elif pathlib.Path(args.output).suffix != ".root":
        print("Output must be a .root file.")
    else:
        main(input_files, args.output)
