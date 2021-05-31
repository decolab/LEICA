function plot_nodes_in_cortex (cortex, V, coord, ori, a, thresh, map, cind, strcont, str, redux)
% 

% Set defaults
if isempty(thresh)
	thresh = 0;	% 0.2;
else
	thresh = thresh/max(abs(V));
end
if isempty(str)
	str = cell(1);
end

% PLOT CORTEX (requires SPM)
cortex.pial = mapPial(cortex.file);
sregion = smooth3(cortex.pial);
psregion = patch(isosurface(sregion,cortex.val,'verbose'), 'FaceColor', cortex.color, 'EdgeColor', 'none');
reducepatch(psregion,redux,'verbose');
isonormals(sregion,psregion);
set(psregion,'FaceAlpha', cortex.transparency);			%transparency

% PLOT NODES
V = V/max(abs(V));
[x,y,z] = sphere;
x = a*x;
y = a*y;
z = a*z;

% Find which nodes have significant changes in in-strength and out-strength
[r, c] = find(str{:,:});

% Plot nodes
for n = 1:length(V)
    if V(n)>0
		if V(n)>thresh && ismember(n,r)
			surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',[1 V(n) 0],'EdgeColor',sum(cind(c(r==n),:),1),'FaceAlpha',0.5);
		elseif V(n)>thresh
			surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',[1 V(n) 0],'EdgeColor','none','FaceAlpha',0.5);
        else
			surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',[1 V(n) 0],'EdgeColor','none','FaceAlpha',0.1);
		end
    elseif V(n)<0
		if abs(V(n))>thresh && ismember(n,r)
			surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',[0 -V(n) 1],'EdgeColor',sum(cind(c(r==n),:),1),'FaceAlpha',0.5);
		elseif abs(V(n))>thresh
			surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',[0 -V(n) 1],'EdgeColor','none','FaceAlpha',0.5);
		else
			surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',[0 -V(n) 1],'EdgeColor','none','FaceAlpha',0.1);
		end
    elseif V(n)==0
        surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',[1 1 1],'EdgeColor','none','FaceAlpha',0.1);
    end
end

% Plot significant connections
n_strong = find(V > thresh); n_weak = find(V < -thresh);
if ~isempty(map) && iscell(map)
	[thr,con] = size(map);
	for r = 1:thr
		for c = 1:con
			for f = 1:numel(map{r,c})
				a(f) = sum(map{r,c}{f},'all');
			end
			a = a./max(a,[],'all','omitnan');
			[~,ai(r,c)] = max(a);
			
			for f = 1:numel(map{r,c})
				[l(:,1), l(:,2)] = find(map{r,c}{f});
				for s = 1:size(l, 1)
					p2 = [coord(l(s,1),2)+ori(1) coord(l(s,1),1)+ori(2) coord(l(s,1),3)+ori(3)];
					p1 = [coord(l(s,2),2)+ori(1) coord(l(s,2),1)+ori(2) coord(l(s,2),3)+ori(3)];

					h(f,c) = mArrow3(p1,p2, 'color',cind(c,:)+[0 1-a(f) 0], 'stemWidth',0.15, 'facealpha',0.7); hold on;
                    % mArrow3(p1,p2, 'color',cind{c}, 'stemWidth',0.15, 'facealpha', a(f));
				end
				clear l
			end
		end
	end
	if ~isempty(strcont)
		for c = 1:con
			l(c) = h(ai(c),c);
		end
		legend(l, strcont, 'Location','southoutside');
	end
elseif ~isempty(map)
	[l(:,1), l(:,2)] = find(map);
	for s = 1:size(l, 1)
		p2 = [coord(l(s,1),2)+ori(1) coord(l(s,1),1)+ori(2) coord(l(s,1),3)+ori(3)];
		p1 = [coord(l(s,2),2)+ori(1) coord(l(s,2),1)+ori(2) coord(l(s,2),3)+ori(3)];
        
		mArrow3(p1,p2, 'color','r', 'stemWidth',0.15, 'facealpha',0.7); hold on;
	end
elseif numel(n_strong)>1 && numel(n_weak)>1
    for a = 1:numel(n_strong)
        n = n_strong(a);
        for b=1:a
            p = n_strong(b);
            c1 = [coord(n,2)+ori(1) coord(n,1)+ori(2) coord(n,3)+ori(3)];
            c2 = [coord(p,2)+ori(1) coord(p,1)+ori(2) coord(p,3)+ori(3)];
            plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)], 'Color','r');	%cmap(IDX(t),:));
        end
    end
    
    for a = 1:numel(n_weak)
        n = n_weak(a);
        for b=1:a
            p = n_weak(b);
            c1 = [coord(n,2)+ori(1) coord(n,1)+ori(2) coord(n,3)+ori(3)];
            c2 = [coord(p,2)+ori(1) coord(p,1)+ori(2) coord(p,3)+ori(3)];
            plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)], 'Color','b');	%cmap(IDX(t),:));
        end
    end
elseif numel(n_strong)>1
    for a = 1:numel(n_strong)
        n = n_strong(a);
        for b=1:a
            p = n_strong(b);
            c1 = [coord(n,2)+ori(1) coord(n,1)+ori(2) coord(n,3)+ori(3)];
            c2 = [coord(p,2)+ori(1) coord(p,1)+ori(2) coord(p,3)+ori(3)];
            plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)], 'Color','r');	%cmap(IDX(t),:));
        end
    end
end

axis off;
axis equal


% -------------------------------------------------------
% Setting image properties - light, material, angle
% -------------------------------------------------------
set(gcf,'Renderer', 'OpenGL') % USE UNDER LINUX FOR TRANSPARENCY 
view(3); axis off;
daspect([1 1 1]);
pbaspect([1 1 1]);
%daspect([0.90  1.37  0.90])
set(gca,'CameraViewAngle', 6);
set(gca, 'Projection', 'orthographic')
set(gca, 'CameraTarget', [51 68 90])
%set(gca, 'CameraPosition', [51 1360 90]) % saggital
%set(gca, 'CameraPosition', [1020 1360 1800]) % free angle
%set(gca, 'CameraUpVector', [1 0 0])

%view([-90 60])
%view([90 -90]) % ventral
view([-90 90]) % top
%view([0 0]) % R sideways
% view([-180 0]) % L sideways
% view([45 20]) % perspective
%view([90 0]) % front

material dull; lighting phong;
camlight;
rotate3d;

end

function pial = mapPial(region)

VG = spm_vol(region(1,:));
pial = zeros(VG.dim(1:3)); 
for i = 1:VG.dim(3)
  pial(:,:,i) = spm_slice_vol(VG,spm_matrix([0 0 i]),VG.dim(1:2),1);
end

end