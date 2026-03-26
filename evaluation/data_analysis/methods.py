# unwanted_methods=["Oriented2-Monge", "OrientedPC-MLS", "OrientedWaveJets", "FO2D", "3DQuadric"]
unwanted_methods=["OrientedPC-MLS", "OrientedWaveJets", "FO2D", "NormCov3D", "VarifoldsMeanPlane", "Oriented2-Monge"]

all_methods = ["PCA", "3DQuadric", "Mean", "Cov2D", "NormCov2D", "NormCov3D", "ShapeOperator", "2-Monge", "Oriented2-Monge","FO2D", "PC-MLS", "OrientedPC-MLS", "JetFitting", "WaveJets", "OrientedWaveJets", "Sphere", "APSS", "UnorientedSphere", "ASO","Varifolds", "VarifoldsMeanPlane", "VCM", "AvgHexagram"]
mean_curvature_estimation=["3DQuadric", "Mean", "Cov2D", "NormCov2D", "NormCov3D", "ShapeOperator", "2-Monge", "Oriented2-Monge","FO2D", "PC-MLS", "OrientedPC-MLS", "JetFitting", "WaveJets", "OrientedWaveJets", "Sphere", "APSS", "UnorientedSphere", "ASO","Varifolds", "VarifoldsMeanPlane", "VCM", "AvgHexagram"]
principal_curvature_estimation=["3DQuadric", "Cov2D", "NormCov2D", "NormCov3D", "ShapeOperator", "2-Monge", "Oriented2-Monge","FO2D", "PC-MLS", "OrientedPC-MLS", "JetFitting", "WaveJets", "OrientedWaveJets", "ASO", "Varifolds", "VarifoldsMeanPlane", "VCM", "AvgHexagram"]
normal_estimation=["PCA", "3DQuadric", "WaveJets", "Mean", "2-Monge", "Oriented2-Monge","FO2D", "PC-MLS", "OrientedPC-MLS", "JetFitting", "Sphere", "APSS", "UnorientedSphere", "ASO", "Varifolds", "VarifoldsMeanPlane"]
curvature_direction=["3DQuadric", "Cov2D", "NormCov2D", "NormCov3D", "ShapeOperator", "2-Monge", "Oriented2-Monge","FO2D", "PC-MLS", "OrientedPC-MLS", "JetFitting", "ASO", "Varifolds", "VarifoldsMeanPlane", "VCM", "AvgHexagram"]

mls_methods=["PC-MLS", "OrientedPC-MLS", "Sphere", "APSS", "UnorientedSphere", "ASO", "FO2D"]

deep_methods=["HoughCNN", "AdaFit", "RefineNet", "DeepFit", "IterNet", "NestiNet", "PCPNet"]

references={
    "PCA" : ("\MtdCovPlane", "Hoffman1987"),
    "Mean": ("\MtdMean", "Pottmann2007"),
    "Cov2D": ("\MtdCovTwoD", "berkmann1994computation,Digne2014"),
    "NormCov2D": ("\MtdNormCovTwoD", "berkmann1994computation,Digne2014"),
    "NormCov3D": ("\MtdNormCovThreeD", "liang1990representation,Digne2014"),
    "NormalW":  ("\MtdNormalW", "Cheng2009"),
    "ShapeOperator": ("\MtdShapeOp", "Kalogerakis2007"),
    "2-Monge": ("\MtdMonge", "Gray1998"),
    # "Oriented2-Monge": ("\MtdMonge", ""),
    "FO2D": ("\MtdPCMLS", ""),
    "PC-MLS": ("\MtdPCMLS", "Ridel2015"),
    "OrientedPC-MLS": ("\MtdPCMLS", ""),
    "JetFitting": ("\MtdJet", "Cazals2005"),
    "WaveJets": ("\MtdWaveJet", "Bearzi2018"),
    "PSS":("\MtdPSS", "Amenta2004"),
    "OrientedWaveJets": ("\MtdWaveJet", ""),
    "Sphere": ("\MtdSphere", "Guennebaud2007"),
    "APSS": ("\MtdAPSS", "Guennebaud2007"),
    "UnorientedSphere": ("\MtdUnorientedSphere", "Chen2013"),
    "ASO": ("\MtdASO", "Lejemble2021"),
    "3DQuadric": ("\MtdQuadricThreeD", "Khameneifar2018"),
    "Varifolds": ("\MtdVarifolds", "Buet2022"),
    "VCM": ("\MtdVCM", "Merigot2011"),
    "AvgHexagram": ("\MtdAvgHex","Lachaud2023"),
    "HoughCNN": ("\\texttt{HoughCNN}", "Boulch2016"),
    "AdaFit": ("\\texttt{AdaFit}", "Zhu2021"),
    # "RefineNet": ("\\texttt{RefineNet}", "Zhou2023RefineNet"),
    "DeepFit": ("\\texttt{DeepFit}", "Yizhak2020"),
    # "IterNet": ("\\texttt{IterNet}", "Boulch2019"),
    # "NestiNet": ("\\texttt{NestiNet}", "Boulch2019"),
    "PCPNet": ("\\texttt{PCPNet}", "Guerrero2018"),
    "NeAF": ("\\texttt{NeAF}", "Li2023NeAF"),
}