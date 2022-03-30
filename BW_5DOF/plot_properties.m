function plot_properties(legends,x_label,y_label,title_fig,index_lt,index_ht,legend_font,fig_font)

if nargin == 6
    legend_font = 16;
    fig_font = 20;
end
if nargin == 7
    fig_font = 20;
end

if isempty(legends)
    set(gca,'LineWidth',2,'FontSize',fig_font,'FontWeight','bold','FontName','Times');
    set(gcf,'Position',[1 1 index_lt*round(2000) index_ht*round(2000)]);
    set(get(gca,'xlabel'),'String',x_label,'FontSize',fig_font','FontWeight','bold','FontName','Times','Interpreter','Latex');
    set(get(gca,'ylabel'),'String',y_label,'FontSize',fig_font','FontWeight','bold','FontName','Times','Interpreter','Latex');
    set(gcf,'color','w'); box on; title(title_fig);
else
    legend(legends,'FontSize',legend_font,'FontWeight','bold')
    set(gca,'LineWidth',2,'FontSize',20,'FontWeight','bold','FontName','Times');
    set(gcf,'Position',[1 1 index_lt*round(2000) index_ht*round(2000)]);
    set(get(gca,'xlabel'),'String',x_label,'FontSize',fig_font','FontWeight','bold','FontName','Times','Interpreter','Latex');
    set(get(gca,'ylabel'),'String',y_label,'FontSize',fig_font','FontWeight','bold','FontName','Times','Interpreter','Latex');
    set(gcf,'color','w'); box on; title(title_fig);
end
end