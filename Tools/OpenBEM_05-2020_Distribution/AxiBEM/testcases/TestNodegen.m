
% This script present several examples of meshes in 2D and Axisymmetrical,
% using the function "nodegen.m". The features of the function are tested.


% Test 1: axisymmetrical sphere testing multiple segments and quad/linear elements
clear
R=1;               % Radius of the sphere
quadelem=1;        % If 0, linear elements (2 nodes) are used instead of quadratic (3 nodes)
segments=[0 R R 0 20 R 0; R 0 0 -R 20 R 0];
[rzb,topology]=nodegen(segments,'y',{},quadelem); % nodes and elements


% Test 2: 2D cylinder testing multiple segments and quad/linear elements
clear
R=1;               % Radius of the sphere
quadelem=1;        % If 0, linear elements (2 nodes) are used instead of quadratic (3 nodes)
segments=[0 R R 0 20 R 0; R 0 0 -R 20 R 0;...
          0 -R -R 0 20 R 0; -R 0 0 R 20 R 0];
[xyb,topology]=nodegen(segments,'y',{},quadelem); % nodes and elements


% Test 3: Ellipse testing segments following inline functions
clear
a=1;b=2;
segmfun{1}=inline('2*sqrt(1-(x/1).^2)');
segmfun{2}=inline('-2*sqrt(1-(x/1).^2)');
segments=[-a 0 a 0 20 j 0; a 0 -a 0 20 2*j 0];
[xyb,topology,rzline,segrzb]=nodegen(segments,'y',segmfun); % nodes and elements


% Test 4: Ellipse testing segments following anonimous functions
clear
a=1;b=2;
segmfun{1}=@(x) b*sqrt(1-(x/a).^2);
segmfun{2}=@(x) -b*sqrt(1-(x/a).^2);
segments=[-a 0 a 0 20 j 0; a 0 -a 0 20 2*j 0];
[xyb,topology,rzline,segrzb]=nodegen(segments,'y',segmfun); % nodes and elements


% Test 5: Ellipse testing segments following anonimous functions, another definition
clear
a=1;b=2;
segmfun{1}=@(y) a*sqrt(1-(y/b).^2);
segmfun{2}=@(y) -a*sqrt(1-(y/b).^2);
segments=[0 b 0 -b 20 -j 0; 0 -b 0 b 20 -2*j 0];
[xyb,topology,rzline,segrzb]=nodegen(segments,'y',segmfun); % nodes and elements


% Test 6: Two rectangles, test of varying length elements
clear
a=4;b=1;
segments=[0 b a b j 0 0; a b a 0 20 0 0; a 0 0 0 j 0 0; 0 0 0 b 3 0 0;...
          0 -b a -b j 0 0; a -b a -2*b 3 0 0; a -2*b 0 -2*b j 0 0; 0 -2*b 0 -b 20 0 0];
[xyb,topology,rzline,segrzb]=nodegen(segments,'y'); % nodes and elements


% Still to be tested: Bézier, varying element length, mesh density,
% multiple bodies

