from typing import Tuple
import tensorflow as tf
from tensorflow.keras import layers

# Cerebras libs
'''
from cerebras_reference_implementations.common.tf.estimator.cs_estimator_spec import (
    CSEstimatorSpec,
)
'''


def unet(input_shape: Tuple[int, int, int]) -> tf.keras.Model:

    input_layer = layers.Input(input_shape)
    x = input_layer

    # Encoder
    x = layers.Conv2D(32, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)
    x = layers.Conv2D(32, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)
    skip1 = x
    x = layers.MaxPool2D(2)(x)

    x = layers.Conv2D(64, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)
    x = layers.Conv2D(64, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)
    skip2 = x
    x = layers.MaxPool2D(2)(x)

    x = layers.Conv2D(128, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)
    x = layers.Conv2D(128, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)
    skip3 = x
    x = layers.MaxPool2D(2)(x)

    # Bottleneck
    x = layers.Conv2D(256, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)
    x = layers.Conv2D(256, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)

    # Decoder
    x = layers.UpSampling2D(2)(x)
    x = layers.Concatenate(axis=-1)([x, skip3])
    x = layers.Conv2D(128, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)
    x = layers.Conv2D(128, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)

    x = layers.UpSampling2D(2)(x)
    x = layers.Concatenate(axis=-1)([x, skip2])
    x = layers.Conv2D(64, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)
    x = layers.Conv2D(64, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)

    x = layers.UpSampling2D(2)(x)
    x = layers.Concatenate(axis=-1)([x, skip1])
    x = layers.Conv2D(32, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)
    x = layers.Conv2D(32, 3, padding='same', kernel_initializer='he_normal')(x)
    x = layers.ReLU()(x)

    x = layers.Conv2D(1, 1, activation='sigmoid', padding='same', kernel_initializer='he_normal')(x)

    model = tf.keras.Model(input_layer, x)
    return model


# Model function
'''
def model_fn(features, labels, mode, params):
    loss, logits = build_model(features, labels, mode, params['model'])

    train_op = None
    host_call = None
    if mode == tf.estimator.ModeKeys.TRAIN:
        optimizer_params = params["optimizer"]

        optimizer_type = optimizer_params.get("optimizer_type", None)
        if optimizer_type is None or optimizer_type.lower() == "adam":
            opt = tf.compat.v1.train.AdamOptimizer(
                learning_rate=optimizer_params['learning_rate'],
                beta1=optimizer_params['beta1'],
                beta2=optimizer_params['beta2'],
                epsilon=optimizer_params['epsilon'],
            )
        elif optimizer_type.lower() == "sgd":
            opt = tf.compat.v1.train.GradientDescentOptimizer(
                learning_rate=optimizer_params["learning_rate"]
            )
        else:
            raise ValueError(f'Unsupported optimizer {optimizer_type}')

        train_op = opt.minimize(
            loss=loss,
            global_step=tf.compat.v1.train.get_or_create_global_step(),
        )
    elif mode == tf.estimator.ModeKeys.EVAL:

        def build_eval_metric_ops(logits, labels):
            return {
                "accuracy": tf.compat.v1.metrics.accuracy(
                    labels=labels, predictions=tf.argmax(input=logits, axis=1),
                ),
            }

        host_call = (build_eval_metric_ops, [logits, labels])
    else:
        raise ValueError("Only TRAIN and EVAL modes supported")

    espec = CSEstimatorSpec(
        mode=mode, loss=loss, train_op=train_op, host_call=host_call,
    )

    return espec
'''


