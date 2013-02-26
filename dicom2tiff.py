import dicom
import Image
meta=dicom.read_file("dicomimage.dcm")
TT=Image.frombuffer("L",imSize,meta.PixelData,"raw","",0,1)
TT.save("testOUTPUT.tiff","TIFF",compression="none")