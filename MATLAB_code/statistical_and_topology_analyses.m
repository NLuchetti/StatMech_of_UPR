%% CONFIGURATION MODELS (cm)

clear all; clc;

sc_mx = readmatrix('matrix.txt'); %Ajacency matrix of native model
for k=1:10
    eval(sprintf('rw%d = sym_generate_srand(sc_mx,1000*%d);',k,k));
    output_name = ['cm_' num2str(k) '.txt'];
    writematrix(eval(sprintf('rw%d',k)), output_name, 'Delimiter', 'tab');
end


%% NETWORK METRICS OF NATIVE MODELS 

clear all; clc;

M1 = readmatrix('matrix.txt');
G1 = graph(M1);

net1 = [degree(G1) centrality(G1, 'closeness') centrality(G1, 'betweenness')];

writematrix(net1,'metrics.txt','Delimiter','tab');


%% NETWORK METRICS OF CONFIGURATION MODELS

clear all; clc;

M = cell(1, 10);
G = cell(1, 10);
net = cell(1, 10);
names = ["degree", "closeness", "betweenness"];
for k=1:10
    filename = ['cm_' num2str(k) '.txt'];
    eval(sprintf('M{%d} = readmatrix(filename);',k));
    G{k} = graph(M{k});
    net{k} = [degree(G{k}) centrality(G{k}, 'closeness') centrality(G{k}, 'betweenness')];
    output_name = ['cm_metrics_' num2str(k) '.txt'];
    writematrix(names,output_name,'Delimiter','tab','WriteMode','append');
    writematrix(net{k},output_name,'Delimiter','tab','WriteMode','append');
end


%% MULTIPLE COMPARISON TEST (example with three models)

clear all; clc;

M1 = readmatrix('matrix1.txt');
M2 = readmatrix('matrix2.txt');
M3 = readmatrix('matrix3.txt');

cl = 'closeness'; bet = 'betweenness';

for k=1:3
    eval(sprintf('G%d = graph(M%d);',k,k));
    eval(sprintf('net%d = [degree(G%d) centrality(G%d, cl) centrality(G%d, bet)];',k,k,k,k));
end
for k=1:3 %Normalized metrics
    eval(sprintf('net%d(:,1) = net%d(:,1)/length(net%d);',k,k,k));
    eval(sprintf('net%d(:,2) = net%d(:,2)*(length(net%d)-1);',k,k,k));
    eval(sprintf('net%d(:,3) = (net%d(:,3)-min(net%d(:,3)))/(max(net%d(:,3))-min(net%d(:,3)));',k,k,k,k,k));
end

clearvars -except net*;

for k=1:3 %Resize arrays
    eval(sprintf('net%d(length(net%d(:,1))+1:lenght_of_longest_array,:)=missing;',k,k));
end


metric1 = [net1(:,1) net2(:,1) net3(:,1)];
metric2 = [net1(:,2) net2(:,2) net3(:,2)];
metric3 = [net1(:,3) net2(:,3) net3(:,3)];

clearvars -except deg cl bet;

[p1,tbl1,stats1] = kruskalwallis(metric1);
[p2,tbl2,stats2] = kruskalwallis(metric2);
[p3,tbl3,stats3] = kruskalwallis(metric3);
[results1,~,~,gnames1] = multcompare(stats1,"CriticalValueType","bonferroni");
[results2,~,~,gnames2] = multcompare(stats2,"CriticalValueType","bonferroni");
[results3,~,~,gnames3] = multcompare(stats3,"CriticalValueType","bonferroni");

clearvars -except result*;


%% GHD evaluation 
clear all; clc;
ghd = cell(1, 10);

M = cell(1,10);
for k = 1:10
    filename = sprintf('cm_%d.txt',k);
    M{k} = load(filename); 
end

sc_mx = readmatrix('matrix.txt');

for l=1:length(M)
    a_M1 = sum(sum(sc_mx));
    a_M2 = sum(sum(sc_mx));
    a1_M1 = [];
    a1_M2 = [];
    h = 1;
    for i=1:length(sc_mx)
        for j=1:length(sc_mx)
            a1_M1(h) = sc_mx(i,j)-a_M1/((length(sc_mx))*(length(sc_mx)-1));
            a1_M2(h) = M{l}(i,j)-a_M2/((length(sc_mx))*(length(sc_mx)-1));
        h=h+1;
        end
    end

    a_11 = [];
    for i=1:length(a1_M1)
        a_11(i) = (a1_M1(i)-a1_M2(i))^2;
    end
    ghd{1,l} = sum(a_11)/((length(sc_mx))*(length(sc_mx)-1));
end

med = [];
for i=1:size(ghd,1)
    med(i) = mean(ghd{i});
end