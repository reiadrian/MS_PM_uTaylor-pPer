function string = f_ProxString_Curva(fId)

   %Encuentra la pr�xima secci�n (busca el pr�ximo string ignorando comentarios y los signo : y =).
   string = '';
   while isempty(string)
      string = textscan(fId,'%s',1,'Delimiter',' :=','MultipleDelimsAsOne',1,'CommentStyle','$');
      string = string{1}{1};
   end