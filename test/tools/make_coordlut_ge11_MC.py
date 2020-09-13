import FWCore.ParameterSet.Config as cms


from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('Whatever',Run3)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')  # load from DB
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')
print "Using GlobalTag: %s" % process.GlobalTag.globaltag.value()

# Fake alignment is/should be ideal geometry
# ==========================================
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")
process.preferFakeAlign = cms.ESPrefer("FakeAlignmentSource")

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

process.analyzer1 = cms.EDAnalyzer("MakeCoordLUTGE11",
    # Verbosity level
    verbosity = cms.untracked.int32(1),

    # Output diectory
    outdir = cms.string("./pc_luts/firmware_ge11_MC/"),

    # Produce "validate.root" to validate the LUTs
    please_validate = cms.bool(True),
)

process.path1 = cms.Path(process.analyzer1)
