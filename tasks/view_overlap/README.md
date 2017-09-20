# View overlap prediction

Learning to predict how much images "overlap" with one another (i.e., form loop closures).

## Data

For each house, there are three main files. The first (xxx_iis.txt) has a line for every pair of images whose visible surfaces overlap.   The second (xxx_vi.txt) has a line for a sampling of vertices on the specified surface mesh indicating which images see it. The third (xxx_iv.txt) has a line for each image indicating a sampling of vertices of the mesh are seen by it.   The other files (xxx_iis_iou.png and (xxx_iis_count.png) have the same information in more convenient matrix forms.  Details of the file formats follow.

xxx_iis.txt

    C matterport_camera_poses/xxx.conf
    II iid1 iid2 iou isect union, count1 count2
    followed by many lines starting with 'II' in the same format

        where II is a keyword (it is capital i twice)
        iid1 and iid2 are zero-based integer indices of images in the conf file,
        count1 is the number of pixels in image1 that see surface points visible in image2,
        count2 is the number of pixels in image2 that see surface points visible in image1,
        isect is min(count1, count2),
        union is count0+count1-isect, and
        iou is isect/union

xxx_vi.txt

    C matterport_camera_poses/xxx.conf
    M poisson_meshes/xxx.ply
    VI vid px py pz nx ny nz n iid1 iid2 ...
    followed by many lines starting with 'VI' in the same format

        where VI is a keyword
        where vid is a zero-based integer index of a vertex in the mesh file
        (px,py,pz) is the position of the vertex,
        (nx,ny,nz) is the normal of the vertex,
        n is the number of images that see it, and
        iid1, iid2, ... is a list of n zero-based image indices that see the vertex

xxx_iv.txt

    C matterport_camera_poses/xxx.conf
    M poisson_meshes/xxx.ply
    IV iid n vid1 vid2 ...
    followed by many lines starting with 'IV' in the same format

        where IV is a keyword
        iid is a zero-based integer index of an image in the conf file
        n is the number of vertices seen by it, and
        vid1, vid2, ... is a list of n zero-based vertex indices seen in the image

xxx_iis_iou.png

    An NxN 16-bit png file

        where N is the number of images in the conf file, and
        pixel (i,j) indicates "intersection over union" (iou) of pixels
        in images i, j multiplied by 1000.

xxx_iis_count.pfm

    An NxN 16-bit pfm file

        where N is the number of images in the conf file, and pixel (i,j) indicates the
        number of pixels in image i that see surface points visible in image j.
        The PFM file format has a header line with "Pf" in ASCII, a second line with
        "<width> <height>" in ASCII, and a third line with "-1.0" in ASCII, followed by
        <width>*<height> 4-byte floats in row-major order (i varies more quickly).

Actually, there are four versions of the image-image overlap files (iis, iip, iig, and iiv),
representing different ways in which overlaps were computed.   In each of the first three files,
overlaps between images i and j are counted by looping over pixels of image i and counting
the number of times the ...

  iis - 3D backprojection of pixel in image i is within 5cm of
        the 3D backprojection of an any pixel in image j.
  iig - 3D backprojection of pixel in image i i falls into the
        same 5-cm-wide grid cell as the backprojection of any pixel in image j.
  iip - projection from pixel in image i to 3D and the back into image j falls on
        a pixel in j with depth +/- 10% of observed depth at that pixel.

In the fourth file (iiv), overlaps between image i and j are counted by looping over
vertices of a mesh and counting the number of vertices visible in both images

  iiv - a vertex of a mesh visible at pixel i is also visible in image j


## Code Requirements
- Install [Torch](http://torch.ch/docs/getting-started.html) on a machine with CUDA GPU
- Install [cuDNN v5](https://developer.nvidia.com/cudnn)
- download pretrained resnet model
    ```
    cd trainingcode
    wget https://d2j0dndfm35trm.cloudfront.net/resnet-50.t7
    ```
- matterport data path: ```$path2matterportdata```
- training testing split:  ```../data/trainval```


## Training
```
cd trainingcode
dataRoot=$path2matterportdata loseType=diff th train.lua 
```
This program trains a view overlap prediction model with following options:
```loseType= diff ```  train with regrassion lose for overlap esitmation.
```loseType= binary ```  train without regrassion lose.
The snapshot models will be saved in
```trainingcode/checkpoints/<$name>/<number_of_iteration>.net```

## Testing
```
name=train1diff  th testFeat.lua
```
This program will load the trained model in folder ```$name```, and compute the image feature for all images in the test scene, the result is saved in ```checkpoints/<$name>/result/....h5```

## Evaluation
The following code will load the computed feature vectors in "testing" step and compute normalized discounted cumulative gain for each testing sequence. 

```
cd trainingcode
matlab &
mpOverlapEval('train1diff')
```
