function LV_RV_assign = determine_LV_RV_by_distance_to_inner_surface(nodes, elems,...
    LV_endo_nodes_list, RV_freewall_nodes_list, LV_dis, RV_dis)

LV_endo_nodes = nodes(LV_endo_nodes_list,2:4);
RV_freewall_nodes = nodes(RV_freewall_nodes_list,2:4);

%% LV = 1;
%% RV = 2

node_assign = ones( [size(nodes, 1) 1]);

for nl = 1 : size(nodes, 1)
    dis_lv_square =  (nodes(nl,2) - LV_endo_nodes(:,1)).^2 + ...
                     (nodes(nl,3) - LV_endo_nodes(:,2)).^2 + ...
                     (nodes(nl,4) - LV_endo_nodes(:,3)).^2;
         
    dis_rv_square =   (nodes(nl,2) - RV_freewall_nodes(:,1)).^2 + ...
                      (nodes(nl,3) - RV_freewall_nodes(:,2)).^2 + ...
                      (nodes(nl,4) - RV_freewall_nodes(:,3)).^2 ;
         
    dis_rv_min = min(dis_rv_square);
    dis_lv_min = min(dis_lv_square);
    
    if dis_rv_min <= RV_dis*RV_dis
        node_assign(nl,1) = 2;
        if dis_rv_min > dis_lv_min %%close to LV then belong to LV
            node_assign(nl,1) = 1;
        end
    elseif dis_lv_min > LV_dis*LV_dis && dis_rv_min < 2*RV_dis*RV_dis
        node_assign(nl,1) = 2;
    end
    
end



elem_assign = ones( [size(elems, 1) 1]);
%% now loop each element and decide what to do 
for el = 1 : size(elems,1)
    el_nodes_list =  elems(el, 2:5);
    node_assign_T = node_assign(el_nodes_list);
    if mean(node_assign_T)>=1.5
        elem_assign(el,1) = 2;
    end

end

LV_RV_assign.node_assign = node_assign;
LV_RV_assign.elem_assign = elem_assign;

























