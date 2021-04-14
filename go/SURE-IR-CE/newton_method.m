            
%%
% 通用牛顿法调用及设置说明：
%                                                                                                                                                                                                                                   
% 牛顿法模块 function_name：目标函数名 
% 调用格式为： newton_method（set，目标函数名'，初始点向量，settings_1，settings_2，settings_3，settings_4，settings_5，settings_6）
%         或  newton_method（set，目标函数名')
% set位数大于等于1 
% set（1）：精度参量 必要输入，作为牛顿法收敛精度的参考
% set（2）：是否在执行牛顿法时调用拟牛顿法 默认为 0
%（在维度比较大的时候，或超维面比较复杂，与超二次曲面相差较大时 拟牛顿法可能会有一定优势，可根据实际运算中观测情况自行调整。）
% set（3）：是否在程序运行中显示运行中的关键参数 
% 其中 若初始点未付值 则默认为上次程序未结束而退出 自动装载上次运行保存的信号
% 各组 settings 要与目标函数一致 ， “目标函数名” 与 “精度” 是必要输入
% 其中牛顿法为了保证遇到特殊状况而造成程序 break
% 会实时保存每次运算的关键参量，保存在newtonmethod_coefficients结构体中
%
% newtonmethod_coefficients结构体 储存参量说明：
% function_name: '目标函数名'
% stepindex:      微分法步长，其维度与初始点相同；
% accuricy:       精度参考值
% coefficients:   [第一位是目标函数此时的结果  后面为相应自变量]
% executed_statement: '根据函数名自动生成的执行语句'
% settings_1，settings_2，settings_3，settings_4，settings_5，settings_6.保存当前目标函数的预设参量
% finished ： 用于判断牛顿法程序是否完成
% 
% alarm: 
% 其中初始点为行向量
% 牛顿法程序内部对目标函数调用格式为：目标函数（初始点，settings_1，settings_2....）
% 其中最多支持6组settings，其中settings是目标函数中除了自变量外其他设置参数， settings不应为自变量；

function [best_point]=newton_method(set,function_name,initial_point,settings_1,settings_2,settings_3,settings_4,settings_5,settings_6)
%% 初始化 & 检错保护 

global cont_newton_iterations
global cont_newton_up
global cont_quasi_iterations
global cont_quasi_up

cont_newton_iterations=0;
cont_newton_up=0;
cont_quasi_iterations=0;
cont_quasi_up=0;

global newtonmethod_search_model
% newtonmethod_search_model=1;
if nargin<2            %缺少输入变量
    error('not enough inputs')
end
if nargin>9
    error('too many input arguments')
end
switch length(set)
    case 1
        newtonmethod_search_model=0;
        newtonmethod_whether_display=0;
    case 2
        newtonmethod_search_model=set(2)
        newtonmethod_whether_display=0;
    case 3
         newtonmethod_search_model=set(2)
         newtonmethod_whether_display=set(3)
    otherwise
         error('unexpected settings')
end
accuricy=set(1);
if isempty(function_name)~=1&&ischar(function_name)
    eval(['target_function=@' num2str(function_name) ';']);
else
    error('unexpected function name')
end

if nargin==2        %上次运算未完成，装载已存函数
    load newtonmethod_coefficients
    if strcmp(newtonmethod_coefficients.function_name,function_name)~=1
        error('this inputed function name is different from the previous one')
    end
elseif nargin>2
    if isempty(initial_point)
        error('please input initial point')
    else
        newtonmethod_coefficients.function_name=function_name;
        newtonmethod_coefficients.stepindex=ones(1,length(initial_point))/20;
        newtonmethod_coefficients.accuricy=accuricy;
        newtonmethod_coefficients.finished=0;
        switch nargin
            case 3
                 newtonmethod_coefficients.coefficients=[0  initial_point];
                 newtonmethod_coefficients.coefficients(1)= target_function(initial_point);
                 newtonmethod_coefficients.executed_statement=['target_function(intermediate_variable)'];
            case 4
                 newtonmethod_coefficients.coefficients=[0  initial_point];
                 newtonmethod_coefficients.coefficients(1)= target_function(initial_point,settings_1);
                 newtonmethod_coefficients.executed_statement=['target_function(intermediate_variable,newtonmethod_coefficients.settings_1)'];
                  newtonmethod_coefficients.settings_1=settings_1;
                 
            case 5
                 newtonmethod_coefficients.coefficients=[0  initial_point];
                 newtonmethod_coefficients.coefficients(1)= target_function(initial_point,settings_1,settings_2);
                 newtonmethod_coefficients.executed_statement=['target_function(intermediate_variable,newtonmethod_coefficients.settings_1,newtonmethod_coefficients.settings_2)'];
                  newtonmethod_coefficients.settings_1=settings_1;
                 newtonmethod_coefficients.settings_2=settings_2;
            case 6
                 newtonmethod_coefficients.coefficients=[0  initial_point];
                 newtonmethod_coefficients.coefficients(1)= target_function(initial_point,settings_1,settings_2,settings_3);
                 newtonmethod_coefficients.executed_statement=['target_function(intermediate_variable,newtonmethod_coefficients.settings_1,newtonmethod_coefficients.settings_2,newtonmethod_coefficients.settings_3)'];
                  newtonmethod_coefficients.settings_1=settings_1;
                 newtonmethod_coefficients.settings_2=settings_2;
                 newtonmethod_coefficients.settings_3=settings_3;
            case 7
                 newtonmethod_coefficients.coefficients=[0  initial_point];
                 newtonmethod_coefficients.coefficients(1)= target_function(initial_point,settings_1,settings_2,settings_3,settings_4);
                 newtonmethod_coefficients.executed_statement=['target_function(intermediate_variable,newtonmethod_coefficients.settings_1,newtonmethod_coefficients.settings_2,newtonmethod_coefficients.settings_3,newtonmethod_coefficients.settings_4)'];     
                 newtonmethod_coefficients.settings_1=settings_1;
                 newtonmethod_coefficients.settings_2=settings_2;
                 newtonmethod_coefficients.settings_3=settings_3;
                 newtonmethod_coefficients.settings_4=settings_4;
            case 8
                 newtonmethod_coefficients.coefficients=[0  initial_point];
                 newtonmethod_coefficients.coefficients(1)= target_function(initial_point,settings_1,settings_2,settings_3,settings_4,settings_5);
                 newtonmethod_coefficients.executed_statement=['target_function(intermediate_variable,newtonmethod_coefficients.settings_1,newtonmethod_coefficients.settings_2,newtonmethod_coefficients.settings_3,newtonmethod_coefficients.settings_4,newtonmethod_coefficients.settings_5)'] ;    
                 newtonmethod_coefficients.settings_1=settings_1;
                 newtonmethod_coefficients.settings_2=settings_2;
                 newtonmethod_coefficients.settings_3=settings_3;
                 newtonmethod_coefficients.settings_4=settings_4;
                 newtonmethod_coefficients.settings_5=settings_5;
            case 9
                 newtonmethod_coefficients.coefficients=[0  initial_point];
                 newtonmethod_coefficients.coefficients(1)= target_function(initial_point,settings_1,settings_2,settings_3,settings_4,settings_5,settings_6);
                 newtonmethod_coefficients.executed_statement=['target_function(intermediate_variable,newtonmethod_coefficients.settings_1,newtonmethod_coefficients.settings_2,newtonmethod_coefficients.settings_3,newtonmethod_coefficients.settings_4,newtonmethod_coefficients.settings_5,newtonmethod_coefficients.settings_6)'] ;    
                 newtonmethod_coefficients.settings_1=settings_1;
                 newtonmethod_coefficients.settings_2=settings_2;
                 newtonmethod_coefficients.settings_3=settings_3;
                 newtonmethod_coefficients.settings_4=settings_4;
                 newtonmethod_coefficients.settings_5=settings_5;
                 newtonmethod_coefficients.settings_6=settings_6; 
        end
        save newtonmethod_coefficients newtonmethod_coefficients
    end
end



%% begin newton method

while 1

    cont_newton_iterations=cont_newton_iterations+1;
    load newtonmethod_coefficients
    equ=newtonmethod_coefficients.executed_statement;
    g=zeros(1,length(newtonmethod_coefficients.coefficients)-1);
    G=zeros((length(newtonmethod_coefficients.coefficients)-1),(length(newtonmethod_coefficients.coefficients)-1));
    initial_point=newtonmethod_coefficients.coefficients(2:end);
    initial_result=newtonmethod_coefficients.coefficients(1);
    N=length(initial_point);
    intermediate_variable=initial_point; 
	result_for_break(1,1)=initial_result;
    result_for_break(2:1+N,1)=initial_point';
    for i_1=1:N
        
        intermediate_variable=initial_point;
        intermediate_variable(i_1)=intermediate_variable(i_1)+ newtonmethod_coefficients.stepindex(i_1);
        result_adi= eval([num2str(equ)]);
        result_for_break(1,end+1)=result_adi;
        result_for_break(2:end,end)=intermediate_variable';
        
        intermediate_variable=initial_point;
        intermediate_variable(i_1)=intermediate_variable(i_1)- newtonmethod_coefficients.stepindex(i_1);
        result_mi= eval([num2str(equ)]); 
        result_for_break(1,end+1)=result_mi;
        result_for_break(2:end,end)=intermediate_variable'
        
        g(i_1) = (result_adi-result_mi)/2/newtonmethod_coefficients.stepindex(i_1);
        G(i_1,i_1) = (result_adi+result_mi-initial_result*2)/(newtonmethod_coefficients.stepindex(i_1)^2);
initial_result-result_adi
initial_result-result_mi

        while (sign(initial_result-result_adi)==sign(initial_result-result_mi))&&( sign(initial_result-result_mi)~=0&&sign(initial_result-result_adi)~=0)
            newtonmethod_coefficients.stepindex(i_1)=newtonmethod_coefficients.stepindex(i_1)/10
            save newtonmethod_coefficients newtonmethod_coefficients

            intermediate_variable=initial_point;
            intermediate_variable(i_1)=intermediate_variable(i_1)+ newtonmethod_coefficients.stepindex(i_1);
            result_adi= eval([num2str(equ)]);
            result_for_break(1,end)=result_adi;
            result_for_break(2:end,end)=intermediate_variable';
            intermediate_variable=initial_point;
            intermediate_variable(i_1)=intermediate_variable(i_1)- newtonmethod_coefficients.stepindex(i_1);
            result_mi= eval([num2str(equ)]); 
            result_for_break(1,end)=result_mi;
            result_for_break(2:end,end)=intermediate_variable';

initial_result-result_adi
initial_result-result_mi
newtonmethod_coefficients.stepindex       
            
            g(i_1) = (result_adi-result_mi)/2/newtonmethod_coefficients.stepindex(i_1);
            G(i_1,i_1) = (result_adi+result_mi-initial_result*2)/(newtonmethod_coefficients.stepindex(i_1)^2);
            if sign(initial_result-result_adi)~=sign(initial_result-result_mi)|| sign(initial_result-result_mi)==0||sign(initial_result-result_adi)==0
                  break
            end

        end

    end
   
    for i_1=1:N-1 %生成 hessian 阵
        for i_2=(i_1+1):N
            
        intermediate_variable=initial_point;
        intermediate_variable(i_1)=intermediate_variable(i_1)+ newtonmethod_coefficients.stepindex(i_1);
        intermediate_variable(i_2)=intermediate_variable(i_2)+ newtonmethod_coefficients.stepindex(i_2);
        result_x= eval([num2str(equ)]);
        result_for_break(1,end+1)=result_x;
        result_for_break(2:end,end)=intermediate_variable';
        
        G(i_1,i_2)=result_x;

        intermediate_variable=initial_point;
        intermediate_variable(i_1)=intermediate_variable(i_1)- newtonmethod_coefficients.stepindex(i_1);
        intermediate_variable(i_2)=intermediate_variable(i_2)- newtonmethod_coefficients.stepindex(i_2);
        result_x = eval([num2str(equ)]);
        result_for_break(1,end+1)=result_x;
        result_for_break(2:end,end)=intermediate_variable';
         G(i_1,i_2)= G(i_1,i_2)+result_x;
        
        
        intermediate_variable=initial_point;
        intermediate_variable(i_1)=intermediate_variable(i_1)+ newtonmethod_coefficients.stepindex(i_1);
        intermediate_variable(i_2)=intermediate_variable(i_2)- newtonmethod_coefficients.stepindex(i_2);
        result_x = eval([num2str(equ)]);
        result_for_break(1,end+1)=result_x;
        result_for_break(2:end,end)=intermediate_variable';
         G(i_1,i_2)= G(i_1,i_2)-result_x;

        intermediate_variable=initial_point;
        intermediate_variable(i_1)=intermediate_variable(i_1)- newtonmethod_coefficients.stepindex(i_1);
        intermediate_variable(i_2)=intermediate_variable(i_2)+ newtonmethod_coefficients.stepindex(i_2);
        result_x = eval([num2str(equ)]);
        result_for_break(1,end+1)=result_x;
        result_for_break(2:end,end)=intermediate_variable'
        
        G(i_1,i_2)= G(i_1,i_2)-result_x;
        G(i_1,i_2) = G(i_1,i_2)/newtonmethod_coefficients.stepindex(i_1)/newtonmethod_coefficients.stepindex(i_2)/4;
        G(i_2,i_1)= G(i_1,i_2);
        end
    end
    newtonstep=-G^(-1)*g';
    newtonstep=newtonstep'
    newtonstep=newtonstep*sqrt(newtonmethod_coefficients.stepindex*newtonmethod_coefficients.stepindex')/(newtonstep*newtonstep');

    zoom_step_index=1; %微分步长放缩标记
    cont_1=1;
    cont_for_break=0;
    max_step_scale=1;
    step_scale=1;
%沿着牛顿步长搜索  

        while 1
            load newtonmethod_coefficients
%             step_scale
%             g
            intermediate_variable=newtonmethod_coefficients.coefficients(2:end) + newtonstep*step_scale;
            result_x = eval([num2str(equ)]);

            if result_x>newtonmethod_coefficients.coefficients(1)
                disp(' change ')
                newtonmethod_coefficients.coefficients(1)=result_x
                newtonmethod_coefficients.coefficients(2:end)=intermediate_variable;
                save newtonmethod_coefficients newtonmethod_coefficients
                cont_newton_up=cont_newton_up+1;
                zoom_step_index=0;
                step_scale = step_scale*(1+cont_1);
                cont_1=cont_1+1;
                cont_for_break=0;
                if abs(step_scale)>max_step_scale
                    max_step_scale=abs(step_scale);
                end
            else
                if cont_1==0
                cont_for_break=cont_for_break+1;
                end
                step_scale = -1*sign(step_scale)*(1+abs(step_scale)/10);
                cont_1=0;
            end
            
            if cont_for_break>1
                break
            end
        
        end

        if zoom_step_index==1
            newtonmethod_coefficients.stepindex=newtonmethod_coefficients.stepindex/100;
            save newtonmethod_coefficients newtonmethod_coefficients
        end    
        
        result_x_index=find(result_for_break(1,1:end)==max(result_for_break(1,1:end)));  % 微分法搜索点大于牛顿法时
        if  result_for_break(1,result_x_index(1,1))>newtonmethod_coefficients.coefficients(1,1)
            newtonmethod_coefficients.coefficients=result_for_break(:,result_x_index(1,1))';
            save newtonmethod_coefficients newtonmethod_coefficients
        end

        if isempty( find( result_for_break(1,:)>newtonmethod_coefficients.coefficients(1)))&&(abs(sum(newtonmethod_coefficients.stepindex))<accuricy)
            best_point = newtonmethod_coefficients.coefficients(2:end);
            newtonmethod_coefficients.finished = 1;
            
            save newtonmethod_coefficients newtonmethod_coefficients
            
            break
        end
        clear result_for_break

    %% quasi newton method
newtonmethod_search_model=0

    if newtonmethod_search_model&&zoom_step_index==0;    
        Hk = G^(-1);
        quasi_divition_step=newtonmethod_coefficients.stepindex;
        cont_quasi_negative=0;
while 1

        cont_quasi_iterations=cont_quasi_iterations+1;
        load newtonmethod_coefficients
        result_for_break(1,1)=newtonmethod_coefficients.coefficients(1,1);
        result_for_break(2:1+N,1)=newtonmethod_coefficients.coefficients(1,2:end)';
        Sk=newtonmethod_coefficients.coefficients(2:end)-initial_point;
       
       
        for i_1=1:N % 拟牛顿法求梯度 %

            intermediate_variable=newtonmethod_coefficients.coefficients(2:end);
            intermediate_variable(i_1)=intermediate_variable(i_1)+ quasi_divition_step(i_1);
            result_adi= eval([num2str(equ)]);
            result_for_break(1,end+1)=result_adi;
            result_for_break(2:end,end)=intermediate_variable'

            intermediate_variable=newtonmethod_coefficients.coefficients(2:end);
            intermediate_variable(i_1)=intermediate_variable(i_1)- quasi_divition_step(i_1);
            result_mi= eval([num2str(equ)]); 
            result_for_break(1,end+1)=result_mi;
            result_for_break(2:end,end)=intermediate_variable'
            g1_quasi(1,i_1) = (result_adi-result_mi)/2/quasi_divition_step(i_1);

            result_for_break(1,1)-result_adi
            result_for_break(1,1)-result_mi
            
            
            while sign( result_for_break(1,1)-result_adi)==sign( result_for_break(1,1)-result_mi)&&(sign( result_for_break(1,1)-result_adi)~=0&&sign( result_for_break(1,1)-result_mi)~=0)
                disp('ee')
                disp('ee')
                quasi_divition_step(i_1)=quasi_divition_step(i_1)/20
                newtonmethod_coefficients.stepindex(i_1) = quasi_divition_step(i_1);
                save newtonmethod_coefficients newtonmethod_coefficients
                intermediate_variable=newtonmethod_coefficients.coefficients(2:end);
                intermediate_variable(i_1)=intermediate_variable(i_1)+ quasi_divition_step(i_1);
                result_adi= eval([num2str(equ)]);

                result_for_break(1,end)=result_adi;
                result_for_break(2:end,end)=intermediate_variable';
                intermediate_variable=newtonmethod_coefficients.coefficients(2:end);
                intermediate_variable(i_1)=intermediate_variable(i_1)- quasi_divition_step(i_1);
                result_mi= eval([num2str(equ)]); 
                result_for_break(1,end)=result_mi;
                result_for_break(2:end,end)=intermediate_variable';
                result_for_break(1,1)-result_adi
                result_for_break(1,1)-result_mi
                g1_quasi(1,i_1) = (result_adi-result_mi)/2/quasi_divition_step(i_1);
                if sign( result_for_break(1,1)-result_adi)~=sign( result_for_break(1,1)-result_mi)|| sign( result_for_break(1,1)-result_adi)==0||sign( result_for_break(1,1)-result_mi)==0
                      break
                end
            end
        end

        Yk=g1_quasi-g;

        Hk_1=Hk+((Sk'-Hk*Yk')*(Sk'-Hk*Yk')')/((Sk'-Hk*Yk')'*Yk');
%       Hk_1=((Sk'-Hk*Yk')*(Sk'-Hk*Yk')')/((Sk'-Hk*Yk')'*Yk');
%       Hk_1=Hk+(Sk'*Sk)/(Sk*Sk')-(Hk*Yk'*Yk*Hk)/(Yk*Hk*Yk')
        quasi_step=-Hk_1*g1_quasi';
        quasi_step=quasi_step'
        quasi_step=quasi_step*sqrt((newtonmethod_coefficients.stepindex*newtonmethod_coefficients.stepindex')/(quasi_step*quasi_step'))

        zoom_step_index=1; %微分步长放缩标记
        cont_1=1;
        cont_for_break=0;
        step_scale=1;

          while 1
        disp('inquasi')  
        intermediate_variable=newtonmethod_coefficients.coefficients(2:end)+quasi_step*step_scale;
        result_x = eval([num2str(equ)]);

           if result_x>newtonmethod_coefficients.coefficients(1)
                disp(' change ')
                newtonmethod_coefficients.coefficients(1)=result_x
                newtonmethod_coefficients.coefficients(2:end)=intermediate_variable;
                save newtonmethod_coefficients newtonmethod_coefficients
                cont_quasi_up=cont_quasi_up+1;
                zoom_step_index=0;
                step_scale = step_scale*(1+cont_1);
                cont_1=cont_1+1;
                cont_for_break=0;
                cont_quasi_negative=0;
                 if abs(step_scale)>max_step_scale
                    max_step_scale=abs(step_scale);
                end
            else
                if cont_1==0
                cont_for_break=cont_for_break+1;
                end
                step_scale = -1*sign(step_scale)*(1+abs(step_scale)/10);
                cont_1=0;
           end
            if cont_for_break>1
                break
            end
          end
        if zoom_step_index==1
           cont_quasi_negative=cont_quasi_negative+1;
           quasi_divition_step=quasi_divition_step/10;
           if  isempty( find( result_for_break(1,1:end)>newtonmethod_coefficients.coefficients(1)))
               newtonmethod_coefficients.stepindex=quasi_divition_step*2;
               save newtonmethod_coefficients  newtonmethod_coefficients
           end
        end 
        result_x_index=find(result_for_break(1,1:end)==max(result_for_break(1,1:end)));  % 微分法搜索点大于牛顿法时
        max_step_scale=1;
        if  result_for_break(1,result_x_index(1,1))>newtonmethod_coefficients.coefficients(1,1)
            newtonmethod_coefficients.coefficients=result_for_break(:,result_x_index(1,1))';
            save newtonmethod_coefficients newtonmethod_coefficients
        end
        if  newtonmethod_coefficients.coefficients(1)~=  result_for_break(1,1) %如果最佳点变更，为下一步拟牛顿法初始化
            g=g1_quasi;
            initial_point=newtonmethod_coefficients.coefficients(2:end);
            Hk=Hk_1;
        end
        if ( isempty( find( result_for_break(1,1:end)>newtonmethod_coefficients.coefficients(1)))&&(abs(sum(quasi_divition_step))<accuricy))||cont_quasi_negative>1||sum(isnan(quasi_step))>0
            clear result_for_break
            break
        end

        clear result_for_break
end
   end  
    

    

end
   


end

