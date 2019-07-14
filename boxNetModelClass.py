import sys, os
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
from skimage.transform import resize
from skimage.io import imread
from skimage.util import view_as_windows, pad
import mrcfile

'''
saved_model_cli show  --dir BoxNet2Mask_20180918 --tag_set serve  --signature_def serving_default
'''

class BoxNetPredictor(object):
  '''
  '''
  DESIRED_SAMPLING_RATE=8. # as done in the paper
  MODEL_IMG_SIZE= 256 # as done in the paper
  STRIDES=64 # as done in the paper
  def __init__(self, modelPath="BoxNet2Mask_20180918"):
    '''
    modelPath is a path to the tensorflow
    '''
    #load tf model
    self.sess= tf.Session()
    tf.saved_model.loader.load(self.sess, ["serve"], modelPath)
    graph = tf.get_default_graph()

    self.input_ = graph.get_tensor_by_name('images_predict:0')
    self.output = graph.get_tensor_by_name('softmax_tensor:0')

  def padToRegularSize(self, inputMic, windowSide, stride, fillWith0=True ):
    '''
    pad one mic so that windowing works for all tiles
    '''
    height, width= inputMic.shape[:2]
    paddingHeight= (0, stride- height%stride )
    paddingWidth=  (0, stride- width%stride  )
    
    paddingValues= [paddingHeight, paddingWidth]
    if fillWith0:
      paddedMic= pad(inputMic, paddingValues, mode="constant", constant_values= np.min(inputMic) )
    else:
      paddedMic= pad(inputMic, paddingValues, mode="wrap" )
    return paddedMic, paddingValues

  def _predict(self, oneMic):

    strideSize= BoxNetPredictor.STRIDES
    imgSize= BoxNetPredictor.MODEL_IMG_SIZE
    
    oneMic, paddingTuples= self.padToRegularSize(oneMic, imgSize, strideSize, fillWith0=False)
    windows = view_as_windows(oneMic, (imgSize, imgSize), step=strideSize) #Returns and HxWximgSizeximgSize array
    windowsOriShape = windows.shape #save original dimensions to recover them
    windows = windows.reshape((-1, imgSize, imgSize, 1)) #reshape it to be tensorflow compatible batchSizeximgSizeximgSizex1
    predsList=[]
    for window in windows: #TODO: This loop is inefficient as each tile is processes independently
      window= np.expand_dims( window, 0)
      output_pred= self.sess.run(self.output, feed_dict={self.input_: window})
      output_pred= output_pred.reshape(-1, BoxNetPredictor.MODEL_IMG_SIZE, BoxNetPredictor.MODEL_IMG_SIZE, 3)
      predsList.append( np.squeeze(output_pred))
    predsList= np.stack(predsList, axis=0)

    predsList= predsList.reshape(windowsOriShape+(3,) )
    output= np.ones(oneMic.shape+(3,))*1e-5
    weights= np.zeros_like(output)
    
    for i in range(predsList.shape[0]):
      actual_i= i*strideSize
      for j in range(predsList.shape[1]):
        actual_j= j*strideSize
        output[actual_i:actual_i+imgSize, actual_j:actual_j+imgSize]+= predsList[i,j]
        weights[actual_i:actual_i+imgSize, actual_j:actual_j+imgSize]+= 1
    output[...,0]=0
    output= output/weights
    output = output[paddingTuples[0][0]:-paddingTuples[0][1], paddingTuples[1][0]:-paddingTuples[1][1], ...]
    return output

  def getDownFactor(self, samplingRate):
    return float(samplingRate)/BoxNetPredictor.DESIRED_SAMPLING_RATE
    
  def loadImg(self, fname):
    if fname.endswith(".mrc"):
      with mrcfile.open(fname, permissive=True) as f:
        x= f.data.copy()
    else:
      x= imread(fname)
    x= x.astype(np.float32)
    x= (x-x.mean())/(x.std())
    return x
    
  def predict(self, fname, samplingRate):
    image= self.loadImg(fname)
    downFact= self.getDownFactor(samplingRate)
    x= resize(image, [int(s*downFact) for s in image.shape if s>1])
    mask= self._predict(x)
    return resize(mask, image.shape), image
    
  def close(self):
    self.sess.close()

def test_executeOneMic():
  predictor= BoxNetPredictor()
  micPath= os.path.expanduser( sys.argv[1])
  samplingRate= float( sys.argv[2] )
  pred, img= predictor.predict(micPath, samplingRate)
  fig= plt.figure()
  fig.add_subplot(121)
  plt.imshow(img, cmap="gray")
  fig.add_subplot(122)
  plt.imshow(pred)
  plt.show()

def executeOneProject(micsPath, masksOutPath):
  predictor= BoxNetPredictor()
  for fname in os.listdir(micsPath):
    if not fname.endswith(".mrc"): continue
    pred, img= predictor.predict(os.path.join(micsPath, fname), samplingRate)
#    pred= pred[...,-1] #this is to select only contamination prediction
    print(pred.shape, img.shape)
    fnameOut= os.path.join(masksOutPath, fname)
    print(fnameOut)
    with mrcfile.new(fnameOut) as f:
      f.set_data(pred.astype(np.float32))
  BoxNetPredictor.close()
      
if __name__=="__main__":
  samplingRate= float( sys.argv[1] )
  micsPath= os.path.expanduser( sys.argv[2])
  masksOutPath= os.path.expanduser( sys.argv[3])
  executeOneProject(micsPath, masksOutPath)
  
