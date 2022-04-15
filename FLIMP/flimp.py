
# +
# Importing Libraries
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import backend as K
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv2D, MaxPool2D, Dense, Flatten, Dropout, BatchNormalization
import random as python_random
import tarfile
import urllib
import os
import numpy as np
from skimage.transform import resize
from tensorflow.keras import regularizers

#from matplotlib import pyplot as plt


print('TensorFlow version:', tf.__version__)
print('Is using GPU?', tf.config.list_physical_devices('GPU'))


data_GT = np.genfromtxt('C:/Users/zrw44597/Work/FLImP/data_train_combined.csv', dtype=None, delimiter=',',names=True, unpack = True, encoding=None) 
ground_truth_train = data_GT['z_difference']
image_name_train = data_GT['Image_name']

data_GT_valid = np.genfromtxt('C:/Users/zrw44597/Work/FLImP/data_valid_combined.csv', dtype=None, delimiter=',',names=True, unpack = True, encoding=None)
ground_truth_test =data_GT_valid['z_difference'] 
image_name_test = data_GT_valid['Image_name']

print(len(ground_truth_train))
# +
## Get random batch
image_dir = 'C:/Users/zrw44597/Work/FLImP/Input_diff_all'

def get_random_batch(ground_truth, image_name, batch_size=4):
    #all_keys = list(annot.keys())
    total_examples = len(ground_truth)
    indices = np.random.choice(range(total_examples), batch_size)
    #x = np.zeros((batch_size, 684, 428,3))
    x = np.zeros((batch_size,224,224,1))
    y = np.zeros((batch_size, 1))
    
    for i, index in enumerate(indices):
        image = tf.keras.preprocessing.image.load_img(os.path.join(image_dir, image_name[index]),color_mode="grayscale")
        
        #print(np.shape(image))
        image = tf.image.random_crop(value=np.array(image), size=(224,224))
        arr = tf.keras.preprocessing.image.img_to_array(image)
        #arr = tf.keras.preprocessing.image.smart_resize(arr,[299, 299])
        arr = tf.keras.applications.mobilenet_v2.preprocess_input(arr)
        arr = np.expand_dims(arr, axis=0)
        x[i] = arr
        y[i] = ground_truth[index]
    
    return x, y


def create_data_generator(batch_size, ground_truth, image_name):
    def data_generator():
        while True:
            x, y = get_random_batch(ground_truth, image_name, batch_size)
            yield (x, y)
    return data_generator
# +
batch_size= 16*4
data_generator_train = create_data_generator(batch_size, ground_truth_train, image_name_train)
train_dataset = tf.data.Dataset.from_generator(data_generator_train,output_types=(tf.float64, tf.float64),
                                                output_shapes=(tf.TensorShape([None, None, None, None]),
                                                        tf.TensorShape([None, None])))

data_generator_valid = create_data_generator(batch_size, ground_truth_test, image_name_test)
valid_dataset = tf.data.Dataset.from_generator(data_generator_valid,output_types=(tf.float64, tf.float64),
                                                output_shapes=(tf.TensorShape([None, None, None, None]),
                                                        tf.TensorShape([None, None])))

# +
strategy = tf.distribute.MirroredStrategy()

train_dist_dataset = strategy.experimental_distribute_dataset(train_dataset)
valid_dist_dataset = strategy.experimental_distribute_dataset(valid_dataset)


def my_rmse(y_true, y_pred):
    error = y_true-y_pred    
    sqr_error = K.square(error)
    mean_sqr_error = K.mean(sqr_error)
    sqrt_mean_sqr_error = K.sqrt(mean_sqr_error)
    return sqrt_mean_sqr_error

def resblock(x, kernelsize, filters):
    fx = layers.Conv2D(filters, kernelsize, activation='relu', padding='same')(x)
    #fx = layers.BatchNormalization()(fx)
    fx = layers.Conv2D(filters, kernelsize, padding='same')(fx)
    out = layers.Add()([x,fx])
    out = layers.ReLU()(out)
    #out = layers.BatchNormalization()(out)
    return out

# optimizer, loss, metrics
adam = tf.keras.optimizers.Adam(learning_rate=0.001)

with strategy.scope():

    input_size = (224,224, 1)
    
    image_input = layers.Input(shape=input_size)
    x = layers.Conv2D(4, kernel_size=(3, 3), activation='relu', padding='same', data_format='channels_last')(image_input)
    x = layers.BatchNormalization()(x)
    x = layers.MaxPool2D(pool_size=(2, 2))(x)

    x = layers.Conv2D(8, kernel_size=(3, 3), activation='relu', padding='same', data_format='channels_last')(x)
    x = layers.BatchNormalization()(x)
    x = layers.MaxPool2D(pool_size=(2, 2))(x)
    
    x = layers.Conv2D(16, kernel_size=(3, 3), activation='relu', padding='same', data_format='channels_last')(x)
    x = layers.BatchNormalization()(x)
    x = layers.MaxPool2D(pool_size=(2, 2))(x)

    x = layers.Conv2D(32, kernel_size=(3, 3), activation='relu', padding='same', data_format='channels_last')(x)
    x = layers.BatchNormalization()(x)
    x = layers.MaxPool2D(pool_size=(2, 2))(x)

    x = layers.Conv2D(64, kernel_size=(3, 3), activation='relu', padding='same', data_format='channels_last')(x)
    x = layers.BatchNormalization()(x)
    x = layers.MaxPool2D(pool_size=(2, 2))(x)

    x = resblock(x, kernelsize=3, filters=64)
    x = resblock(x, kernelsize=3, filters=64)

    x = Flatten()(x)
    x = Dense(128, activation='relu')(x)
    outputs = Dense(1)(x)

    model = keras.Model(image_input,outputs)

    
    model.summary()
    
    model.compile(optimizer=adam,
              loss=  my_rmse,
              metrics=['MeanAbsoluteError'])

steps_per_epoch = int(len(ground_truth_train)/batch_size)
validation_steps = int(len(ground_truth_test)/batch_size)

print('Steps per epoch:', steps_per_epoch)
print('Validation steps:', validation_steps)

cp_callback = tf.keras.callbacks.ModelCheckpoint("best_model-cnnRes-parallel-all_LRpt001_{epoch:04d}.h5", monitor='loss', verbose=1,
                                                 save_best_only=True, mode='auto',save_freq=steps_per_epoch*5)

H = model.fit(train_dist_dataset, validation_data=valid_dist_dataset,steps_per_epoch=steps_per_epoch, validation_steps=validation_steps, epochs=15, verbose=2,callbacks=[cp_callback])

model.save('model-cnnRes-parallel-all_LRpt001.h5')
np.save('history-cnnRes-parallel-all_LRpt001.npy',H.history)

#x, y, images = get_random_batch(ground_truth_train, image_name_train, batch_size=8)
#x, y = get_random_batch(ground_truth_test, image_name_test, batch_size=8)
#print(np.shape(x))
#preds = model.predict(x)
#print(preds)
#print(y)

