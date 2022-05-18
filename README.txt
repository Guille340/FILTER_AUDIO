Filter Audio Files
===================

MATLAB code for filtering large audio files using the overlap-save method. 
This approach is much faster than MATLAB's FILTER, since it doesn't need to
load the entire file into RAM.

Improvements have been applied to ensure that there is no jump in value 
in the connecting samples between consecutive processed blocks (14 Aug 2021).

This code uses the "Audio Read" package.

[Guillermo Jiménez-Arranz, 14 Aug 2021]






