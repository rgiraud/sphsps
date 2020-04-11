

function spherical_3D_view(img,B,type,varargin)

img = flipud(img);
B = flipud(B);

[h,w,z] = size(img);
r = zeros(h,w);
g = r;
b = r;

label = -2;
if (nargin >3)
    label = varargin{1};
end
color_label = [1 0 0];
if (nargin >4)
    color_label = varargin{2};
end

color_B = 'k';
if strcmp(type,'img')==1
   if (nargin == 4)
      color_B = varargin{1}; 
   end
end

if (label ~= -2)
    th_lab = 6;
else
    th_lab = 1.25;
end

val_0 = 0.975;
switch type
    case 'labels'
        S = img;
        N = max(S(:));
        %If -1 labels
        sp_pos = (S == 0);
        r(sp_pos) = val_0;
        g(sp_pos) = val_0;
        b(sp_pos) = val_0;
        for i=1:N
            sp_pos = (S == i);
                r(sp_pos) = rand(1)/th_lab + 1 - 1/th_lab - 0.01;
                g(sp_pos) = rand(1)/th_lab + 1 - 1/th_lab - 0.01;
                b(sp_pos) = rand(1)/th_lab + 1 - 1/th_lab - 0.01;
        end
        %specific label
        sp_pos = (S == label);
        r(sp_pos) = color_label(1);
        g(sp_pos) = color_label(2);
        b(sp_pos) = color_label(3);
 
        BB = B;
        BB(1,:) = 1;
        BB(end,:) = 1;
        ppbb = double(cat(3,r,g,b).*repmat(~BB,[1 1 z]));
        pp = double(cat(3,r,g,b).*repmat(~B,[1 1 z]));
    case 'img'
        switch color_B
            case 'k'
                pp = img.*repmat(~B,[1 1 3]);
            case 'w'
                pp = min(img+repmat(B,[1 1 3]),1);
            case 'r'
                pp = max(min(img+double(cat(3,B,-B,-B)),1),0);
            case 'g'
                pp = max(min(img+double(cat(3,-B,B,-B)),1),0);
            case 'b'
                pp = max(min(img+double(cat(3,-B,-B,B)),1),0);
        end
end
    

[x,y,z] = sphere(1000);

figure,
surface(x,y,z, 'FaceColor','texturemap', 'EdgeColor','none', 'Cdata', pp); 
view(3)
axis equal
grid on
xlabel('x');
ylabel('y');
zlabel('z');
view(90,40)






