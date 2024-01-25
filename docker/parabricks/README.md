# Parabricks

Patching NVIDIA latest from:

```
nvcr.io/nvidia/clara/clara-parabricks:4.2.1-1
```

We only add bgzip and tabix to the image, updated image is:

```
zinno/parabricks:4.2.1-1b
```

TensorRT version the model was compiled with prevented use on 4.2.1-1, so we are using 4.1.2-1.ultimaoct2

```
nvcr.io/ea-nvidia-clara-parabricks/clara-parabricks:4.1.2-1.ultimaoct2
```
