function plot_nodes_in_cortex (cortex, V, coord, ori, a, thresh, map, cind, lgnd, str, redux)
%   PLOT_NODES_IN_CORTEX represents the network of interest in cortical
% space.
% INPUTS
%   cortex:	structure which defines cortical space
%   V:      node weights which define network(s) of interest.  May be
%           binary or weighted.
%   coord:  coordinates for nodes in brain space
%   ori:    origin of coordinate space
%   nlgnd:  node labels, arranged as Nx1 string array
%   a:      scaling parameter: controls size of node markers in cortex
%   thresh: threshold for highlighting nodes of interest.  May be left
%           empty if no node(s) should be highlighted.  Absolute value:
%           applies to both positive and negative weights.  Any scaling
%           should be done prior to input into function.
%   map:    edges to plot between nodes.  May be left empty if no edges of
%           interest exist.
%   cind:	structure which defines color index.  cind.node denotes node
%           color key; cind.conn denotes edge color key.
%   lgnd:	strings for figure legend
%   str:    array of nodes to highlight (denotes strength changes).  May be
%           left empty.
%   redux:  proportion of surface faces to keep (<1)
% OUTPUTS: plot of network of interest in cortical space.

% Set defaults
if isempty(thresh)
	thresh = 0;	% 0.2;
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
[x,y,z] = sphere;
x = a*x;
y = a*y;
z = a*z;

% Find which nodes have significant changes in in-strength and out-strength
[r, c] = find(str{:,:});

% Locate significant connections
n_strong = find(V > thresh);
n_weak = find(V < -thresh);

% Plot nodes
for n = 1:length(V)
    if V(n)>0
		if V(n)>thresh && ismember(n,r)
			surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',cind.node(1,:), 'EdgeColor',sum(cind.node(c(r==n),:),1),'FaceAlpha',0.5);
		elseif V(n)>thresh
			surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',cind.node(1,:), 'EdgeColor','none','FaceAlpha',0.5);
        else
			surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',cind.node(1,:), 'EdgeColor','none','FaceAlpha',0.1);
		end
    elseif V(n)<0
		if abs(V(n))>thresh && ismember(n,r)
			surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',cind.node(2,:), 'EdgeColor',sum(cind.node(c(r==n),:),1),'FaceAlpha',0.5);
		elseif abs(V(n))>thresh
			surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',cind.node(2,:), 'EdgeColor','none','FaceAlpha',0.5);
		else
			surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',cind.node(2,:), 'EdgeColor','none','FaceAlpha',0.1);
		end
    elseif V(n)==0 || isnan(V(n))
        surf(x+coord(n,2)+ori(1), y+coord(n,1)+ori(2), z+coord(n,3)+ori(3),'FaceColor',[0 0 0],'EdgeColor','none','FaceAlpha',0.1);
    end
end

% Plot significant connections
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
					h(f,c) = mArrow3(p1,p2, 'color',cind.conn(c,:), 'stemWidth',0.15, 'tipWidth',0.8, 'facealpha',0.8); hold on;   % sets arrow color to contrast
				end
				clear l
			end
		end
	end
	if ~isempty(lgnd)
		for c = 1:con
			l(c) = h(ai(c),c);
		end
		legend(l, lgnd, 'Location','southoutside');
	end
elseif ~isempty(map)
	[l(:,1), l(:,2)] = find(map);
	for s = 1:size(l, 1)
		p2 = [coord(l(s,1),2)+ori(1) coord(l(s,1),1)+ori(2) coord(l(s,1),3)+ori(3)];
		p1 = [coord(l(s,2),2)+ori(1) coord(l(s,2),1)+ori(2) coord(l(s,2),3)+ori(3)];
        
		mArrow3(p1,p2, 'color',cind.conn, 'stemWidth',0.15, 'tipWidth',0.8, 'facealpha',0.7); hold on;
	end
elseif numel(n_strong)>1 || numel(n_weak)>1
    for a = 1:numel(n_strong)
        n = n_strong(a);
        for b=1:a
            p = n_strong(b);
            c1 = [coord(n,2)+ori(1) coord(n,1)+ori(2) coord(n,3)+ori(3)];
            c2 = [coord(p,2)+ori(1) coord(p,1)+ori(2) coord(p,3)+ori(3)];
            plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)], 'Color',cind.conn(1,:), 'LineWidth',0.75);
            %cmap(IDX(t),:));
        end
    end
    
    for a = 1:numel(n_weak)
        n = n_weak(a);
        for b=1:a
            p = n_weak(b);
            c1 = [coord(n,2)+ori(1) coord(n,1)+ori(2) coord(n,3)+ori(3)];
            c2 = [coord(p,2)+ori(1) coord(p,1)+ori(2) coord(p,3)+ori(3)];
            plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)], 'Color',cind.conn(2,:), 'LineWidth',0.75);
            %cmap(IDX(t),:));
        end
    end
elseif numel(n_strong)>1
    for a = 1:numel(n_strong)
        n = n_strong(a);
        for b=1:a
            p = n_strong(b);
            c1 = [coord(n,2)+ori(1) coord(n,1)+ori(2) coord(n,3)+ori(3)];
            c2 = [coord(p,2)+ori(1) coord(p,1)+ori(2) coord(p,3)+ori(3)];
            plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)], 'Color','r', 'LineWidth',0.75);
            %cmap(IDX(t),:));
        end
    end
elseif numel(n_weak)>1
    for a = 1:numel(n_weak)
        n = n_weak(a);
        for b=1:a
            p = n_weak(b);
            c1 = [coord(n,2)+ori(1) coord(n,1)+ori(2) coord(n,3)+ori(3)];
            c2 = [coord(p,2)+ori(1) coord(p,1)+ori(2) coord(p,3)+ori(3)];
            plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)], 'Color',cind.conn(2,:), 'LineWidth',0.75);
            %cmap(IDX(t),:));
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