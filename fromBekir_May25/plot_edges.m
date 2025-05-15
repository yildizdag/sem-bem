function  [] = plot_edges(group_ind,posn,elementpoints)

ne = size(group_ind,1); % number of all edges

figure(300)
for i = 1:ne
    
    % Plot edge #1 of 4 edges of the ith element
    plot(posn(elementpoints(i,group_ind(i,:)==1|group_ind(i,:)==2|group_ind(i,:)==5),1),...
        posn(elementpoints(i,group_ind(i,:)==1|group_ind(i,:)==2|group_ind(i,:)==5),2),'-')

    % label edge
    txt1 = mean(posn(elementpoints(i,group_ind(i,:)==1|group_ind(i,:)==2|group_ind(i,:)==5),:));
    text(txt1(1),txt1(2),num2str(5),'fontsize',12)

    hold on

    % Plot edge #2 of 4 edges of the ith element
    plot(posn(elementpoints(i,group_ind(i,:)==4|group_ind(i,:)==7|group_ind(i,:)==3),1),...
    posn(elementpoints(i,group_ind(i,:)==4|group_ind(i,:)==7|group_ind(i,:)==3),2),'-')

    % label edge
    txt2 = mean(posn(elementpoints(i,group_ind(i,:)==4|group_ind(i,:)==7|group_ind(i,:)==3),:));
    text(txt2(1),txt2(2),num2str(7),'fontsize',12)

    % Plot edge #3 of 4 edges of the ith element
    plot(posn(elementpoints(i,group_ind(i,:)==1|group_ind(i,:)==8|group_ind(i,:)==4),1),...
        posn(elementpoints(i,group_ind(i,:)==1|group_ind(i,:)==8|group_ind(i,:)==4),2),'-')

    % label edge
    txt3 = mean(posn(elementpoints(i,group_ind(i,:)==1|group_ind(i,:)==8|group_ind(i,:)==4),:));
    text(txt3(1),txt3(2),num2str(8),'fontsize',12) 

    % Plot edge #4 of 4 edges of the ith element    
    plot(posn(elementpoints(i,group_ind(i,:)==2|group_ind(i,:)==6|group_ind(i,:)==3),1),...
        posn(elementpoints(i,group_ind(i,:)==2|group_ind(i,:)==6|group_ind(i,:)==3),2),'-')

    % label edge
    txt4 = mean(posn(elementpoints(i,group_ind(i,:)==2|group_ind(i,:)==6|group_ind(i,:)==3),:));
    text(txt4(1),txt4(2),num2str(6),'fontsize',12)
    
    % Label the element in the center
    txt = mean([txt1;txt2;txt3;txt4]);
    text(txt(1),txt(2),num2str(i),'fontsize',14,'Color','r')

    axis equal
end
end