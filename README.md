Triangular mesh expression of cells.
This code demonstates that two cells collide with each other.
Swelling ratio of each cell is set to be 0.65.
You can change parameter' values in _parameters.h .

"sig" file is original file format to describe the placement of cells.
You see this is almost "off" file but the additional two final lines are different.

Here is the definition of "sig" file.
## sig ##
SIG // header <br>
4 4 0 1// The number of vertices, faces, edges, and cells<br>
1.0 0.0 1.0// x-,y-,and z- components of vertices<br>
. <br>
.<br>
.<br>
3 0 1 2//Components of a face, "3" indicates the number of vertices that composes the face, subsequent digits indicates indices of the vertices.<br>
3 0 2 3 //anti-clockwise is the outward axis of the face<br>
.<br>
.<br>
.<br>
4 0 1 2 3//Components of a cell, "4" indicates the number of faces that composes the cell, subsequent digits indicates indices of the faces. <br>

## How to make and run ##
mkdir output

mkdir data

make

./a.out start_step (usually "./a.out 0")


## Acknowledge ##
Thanks to Hitoshi Saigo
