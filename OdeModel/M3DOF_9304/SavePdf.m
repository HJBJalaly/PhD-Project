function SavePdf(Name)

set(gcf,'Units','centimeters','PaperUnits','centimeters');
Pos=get(gcf,'Position');
set(gcf,'PaperPosition', [0 0 Pos(3) Pos(4)],'PaperSize', [Pos(3) Pos(4)])
% set(gcf,'PaperUnits',Units);
saveas(gcf,[Name,'.pdf'])