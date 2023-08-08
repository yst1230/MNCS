function acc=yst_acc(f_Train,f_Test,train_label,test_label)

Mdl = fitcsvm(f_Train,train_label,'Standardize',true,'KernelFunction', 'linear','KernelScale','auto');
predict_label=predict(Mdl,f_Test);
acc=numel(find(predict_label-test_label==0))/length(test_label);