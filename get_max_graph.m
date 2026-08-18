function [max_graph,graph_fine,smooth]=get_max_graph(spectrum,graph)
smooth=linspace(min(spectrum(:)),max(spectrum(:)),200);

graph_fine = makima(spectrum,graph,smooth);

max_graph=smooth(graph_fine==max(graph_fine(:)));
max_graph=mean(max_graph,'all');
end