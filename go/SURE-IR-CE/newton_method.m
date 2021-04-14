            
%%
% ͨ��ţ�ٷ����ü�����˵����
%                                                                                                                                                                                                                                   
% ţ�ٷ�ģ�� function_name��Ŀ�꺯���� 
% ���ø�ʽΪ�� newton_method��set��Ŀ�꺯����'����ʼ��������settings_1��settings_2��settings_3��settings_4��settings_5��settings_6��
%         ��  newton_method��set��Ŀ�꺯����')
% setλ�����ڵ���1 
% set��1�������Ȳ��� ��Ҫ���룬��Ϊţ�ٷ��������ȵĲο�
% set��2�����Ƿ���ִ��ţ�ٷ�ʱ������ţ�ٷ� Ĭ��Ϊ 0
%����ά�ȱȽϴ��ʱ�򣬻�ά��Ƚϸ��ӣ��볬�����������ϴ�ʱ ��ţ�ٷ����ܻ���һ�����ƣ��ɸ���ʵ�������й۲�������е�������
% set��3�����Ƿ��ڳ�����������ʾ�����еĹؼ����� 
% ���� ����ʼ��δ��ֵ ��Ĭ��Ϊ�ϴγ���δ�������˳� �Զ�װ���ϴ����б�����ź�
% ���� settings Ҫ��Ŀ�꺯��һ�� �� ��Ŀ�꺯������ �� �����ȡ� �Ǳ�Ҫ����
% ����ţ�ٷ�Ϊ�˱�֤��������״������ɳ��� break
% ��ʵʱ����ÿ������Ĺؼ�������������newtonmethod_coefficients�ṹ����
%
% newtonmethod_coefficients�ṹ�� �������˵����
% function_name: 'Ŀ�꺯����'
% stepindex:      ΢�ַ���������ά�����ʼ����ͬ��
% accuricy:       ���Ȳο�ֵ
% coefficients:   [��һλ��Ŀ�꺯����ʱ�Ľ��  ����Ϊ��Ӧ�Ա���]
% executed_statement: '���ݺ������Զ����ɵ�ִ�����'
% settings_1��settings_2��settings_3��settings_4��settings_5��settings_6.���浱ǰĿ�꺯����Ԥ�����
% finished �� �����ж�ţ�ٷ������Ƿ����
% 
% alarm: 
% ���г�ʼ��Ϊ������
% ţ�ٷ������ڲ���Ŀ�꺯�����ø�ʽΪ��Ŀ�꺯������ʼ�㣬settings_1��settings_2....��
% �������֧��6��settings������settings��Ŀ�꺯���г����Ա������������ò����� settings��ӦΪ�Ա�����

function [best_point]=newton_method(set,function_name,initial_point,settings_1,settings_2,settings_3,settings_4,settings_5,settings_6)
%% ��ʼ�� & ����� 

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
if nargin<2            %ȱ���������
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

if nargin==2        %�ϴ�����δ��ɣ�װ���Ѵ溯��
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
   
    for i_1=1:N-1 %���� hessian ��
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

    zoom_step_index=1; %΢�ֲ����������
    cont_1=1;
    cont_for_break=0;
    max_step_scale=1;
    step_scale=1;
%����ţ�ٲ�������  

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
        
        result_x_index=find(result_for_break(1,1:end)==max(result_for_break(1,1:end)));  % ΢�ַ����������ţ�ٷ�ʱ
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
       
       
        for i_1=1:N % ��ţ�ٷ����ݶ� %

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

        zoom_step_index=1; %΢�ֲ����������
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
        result_x_index=find(result_for_break(1,1:end)==max(result_for_break(1,1:end)));  % ΢�ַ����������ţ�ٷ�ʱ
        max_step_scale=1;
        if  result_for_break(1,result_x_index(1,1))>newtonmethod_coefficients.coefficients(1,1)
            newtonmethod_coefficients.coefficients=result_for_break(:,result_x_index(1,1))';
            save newtonmethod_coefficients newtonmethod_coefficients
        end
        if  newtonmethod_coefficients.coefficients(1)~=  result_for_break(1,1) %�����ѵ�����Ϊ��һ����ţ�ٷ���ʼ��
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

