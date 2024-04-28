"""
base class and functions for model training/inference, and visiulations
"""
import os
import pickle
import torch
import numpy as np
import pandas as pd
import torch.nn as nn
import scipy.misc as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image
from scipy.stats import linregress
from keras.models import load_model
from torchvision import transforms, models
from sklearn.metrics import mean_squared_error


def split_df_to2(df, n):
    '''
    args:
        df: pandas dataframe
        n (int): the number of rows for part one
    return:
        df1: dataframe of part1
        df2: dataframe of part2
    '''
    df1 = df.sample(n)
    df2 = df.loc[df.index[~df.index.isin(df1.index)], :]
    return df1, df2


def split_val_train(df, n_val, n_fold):
    '''
    args:
        df: the dataframe for training and validataion for one category
        n_val (int): number of validation data in each fold
    yield fold, index of training data, index of validation data
    '''
    all_idx = df.index
    if n_val*n_fold > all_idx.shape[0]:
        print('not enough data for %s fold cross-validation!' % n_fold)
        return
    val_idx_st, val_idx_ed = 0, n_val
    for fold in range(1, n_fold+1):
        val_idx = all_idx[val_idx_st: val_idx_ed]
        train_idx = all_idx[~all_idx.isin(val_idx)]
        val_idx_st += n_val
        val_idx_ed += n_val
        yield fold, train_idx, val_idx


class EarlyStopping:
    """Early stops the training if validation loss doesn't improve after a
    given patience."""
    def __init__(self, mn_prefix, patience=20, verbose=True, delta=0,
                 min_loss_cutoff=0):
        """
        Args:
            mn_prefix (str): the prefix of the saved model name.
            patience (int): How long to wait after last time validation loss
                improved. Default: 20
            verbose (bool): If True, prints a message for each validation loss
                improvement. Default: True
            delta (float): Minimum change in the monitored quantity to qualify
                as an improvement. Default: 0
            min_loss_cutoff (float): set the minimum loss cutoff if there is
                no validation process during training. For leaf counting
                problem with MSE as the loss, set 0.81 consdering human error
                is 0.5.
        """
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf
        self.delta = delta
        self.mn_prefix = mn_prefix
        self.min_loss_cutoff = -min_loss_cutoff

    def __call__(self, val_loss, model):
        score = -val_loss

        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
        elif score < self.best_score + self.delta:
            self.counter += 1
            print(f'EarlyStopping counter: {self.counter} out of'
                  f' {self.patience}')
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
            self.counter = 0
            if self.best_score > self.min_loss_cutoff:
                self.early_stop = True

    def save_checkpoint(self, val_loss, model):
        '''Saves model when validation loss decrease.'''
        if self.verbose:
            print(f'Validation loss decreased ({self.val_loss_min:.6f} -->'
                  f' {val_loss:.6f}).  Saving model ...')
        torch.save(model.state_dict(), '%s.pt' % self.mn_prefix)
        self.val_loss_min = val_loss


def image_transforms(input_size=224):
    image_transforms_dict = {
        # Train uses data augmentation
        'train':
        transforms.Compose([
            transforms.ColorJitter(),
            transforms.RandomHorizontalFlip(),
            transforms.Resize(size=(input_size, input_size)),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
            ]),
        'valid':
        transforms.Compose([
            transforms.Resize(size=(input_size, input_size)),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
        ]),
    }
    return image_transforms_dict


def train_model_regression(model, dataloaders, criterion, optimizer,
                           model_name_prefix, patience=10, num_epochs=10,
                           is_inception=False):
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model.to(device)
    train_loss_history, valid_loss_history = [], []
    early_stopping = EarlyStopping(model_name_prefix, patience=patience,
                                   verbose=True, min_loss_cutoff=0.8)
    for epoch in range(num_epochs):
        print('Epoch {}/{}'.format(epoch, num_epochs - 1))
        print('-' * 10)
        model.train()
        train_running_loss = 0.0
        for inputs, labels, _ in dataloaders['train']:
            inputs = inputs.to(device)
            labels = labels.to(device)
            optimizer.zero_grad()
            if is_inception:
                outputs, aux_outputs = model(inputs)
                loss1 = criterion(outputs, labels)
                loss2 = criterion(aux_outputs, labels)
                loss = loss1 + 0.4*loss2
            else:
                outputs = model(inputs)
                loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            train_running_loss += loss.item() * inputs.size(0)
        epoch_loss_train = train_running_loss /\
            len(dataloaders['train'].dataset)
        print('Train Loss: {:.4f}'.format(epoch_loss_train))
        train_loss_history.append(epoch_loss_train)

        if 'valid' in dataloaders:
            model.eval()
            valid_running_loss = 0.0
            for inputs, labels, _ in dataloaders['valid']:
                inputs = inputs.to(device)
                labels = labels.to(device)
                outputs = model(inputs)
                loss = criterion(outputs, labels)
                valid_running_loss += loss.item() * inputs.size(0)
            epoch_loss_valid = valid_running_loss /\
                len(dataloaders['valid'].dataset)
            print('Validation Loss: {:.4f}'.format(epoch_loss_valid))
            valid_loss_history.append(epoch_loss_valid)
            early_stopping(epoch_loss_valid, model)
        else:
            early_stopping(epoch_loss_train, model)

        if early_stopping.early_stop:
            print('Early stopping')
            break
    print('Best val loss: {:4f}'.format(early_stopping.val_loss_min))
    model.load_state_dict(torch.load('%s.pt' % model_name_prefix))
    return model, train_loss_history, valid_loss_history


def set_parameter_requires_grad(model, feature_extracting):
    if feature_extracting:
        for param in model.parameters():
            param.requires_grad = False


def initialize_model(model_name='vgg16', num_classes=1, feature_extract=True,
                     use_pretrained=True, inputsize=224):
    model_ft = None
    input_size = 0

    if model_name == "googlenet":
        """ googlenet
        """
        model_ft = models.googlenet(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.fc.in_features
        model_ft.fc = nn.Linear(num_ftrs, num_classes)
        input_size = inputsize

    elif model_name == "resnet152":
        """ Resnet152
        """
        model_ft = models.resnet152(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.fc.in_features
        model_ft.fc = nn.Linear(num_ftrs, num_classes)
        input_size = inputsize

    elif model_name == "vgg16":
        """ VGG16
        """
        model_ft = models.vgg16(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.classifier[6].in_features
        model_ft.classifier[6] = nn.Linear(num_ftrs, num_classes)
        input_size = inputsize

    elif model_name == "resnet18":
        """ Resnet18
        """
        model_ft = models.resnet18(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.fc.in_features
        model_ft.fc = nn.Linear(num_ftrs, num_classes)
        input_size = inputsize

    elif model_name == "alexnet":
        """ Alexnet
        """
        model_ft = models.alexnet(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.classifier[6].in_features
        model_ft.classifier[6] = nn.Linear(num_ftrs, num_classes)
        input_size = inputsize

    elif model_name == "vgg11_bn":
        """ VGG11_bn
        """
        model_ft = models.vgg11_bn(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.classifier[6].in_features
        model_ft.classifier[6] = nn.Linear(num_ftrs, num_classes)
        input_size = inputsize

    elif model_name == "squeezenet":
        """ Squeezenet
        """
        model_ft = models.squeezenet1_0(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        model_ft.classifier[1] = nn.Conv2d(512, num_classes,
                                           kernel_size=(1, 1),
                                           stride=(1, 1))
        model_ft.num_classes = num_classes
        input_size = inputsize

    elif model_name == "densenet":
        """ Densenet
        """
        model_ft = models.densenet121(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.classifier.in_features
        model_ft.classifier = nn.Linear(num_ftrs, num_classes)
        input_size = inputsize

    elif model_name == "inception":
        """ Inception v3
        Be careful, expects (299,299) sized images and has auxiliary output
        """
        model_ft = models.inception_v3(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        # Handle the auxilary net
        num_ftrs = model_ft.AuxLogits.fc.in_features
        model_ft.AuxLogits.fc = nn.Linear(num_ftrs, num_classes)
        # Handle the primary net
        num_ftrs = model_ft.fc.in_features
        model_ft.fc = nn.Linear(num_ftrs, num_classes)
        input_size = inputsize

    else:
        print("Invalid model name, exiting...")
        exit()

    return model_ft, input_size


def Plot(history_pickle_fn, output_fn):
    """
    plot training process
    You can load the dict back using pickle.load(open('*.p', 'rb'))
    """
    hist = pickle.load(open(history_pickle_fn, 'rb'))
    mpl.rcParams['figure.figsize'] = [7.5, 3.25]
    fig, axes = plt.subplots(nrows=1, ncols=2)
    ax1 = axes[0]
    ax1.plot(hist['acc'])
    ax1.plot(hist['val_acc'])
    ax1.set_title('model accuracy')
    ax1.set_ylabel('accuracy')
    ax1.set_xlabel('epoch')
    ax1.set_ylim(0, 1.01)
    ax1.legend(['train', 'validation'], loc='lower right')
    max_acc = max(hist['val_acc'])

    ax2 = axes[1]
    ax2.plot(hist['loss'])
    ax2.plot(hist['val_loss'])
    ax2.set_title('model loss')
    ax2.set_ylabel('loss')
    ax2.set_xlabel('epoch')
    ax2.legend(['train', 'validation'], loc='upper right')
    plt.tight_layout()
    plt.savefig('%s_%s.png' % (max_acc, output_fn))
    plt.clf()


def hyp_image_to_array(hyp_img_dir):
    '''
    convert hyperspectral image dataset to numpy array

    Returns:
        numpy array object with shape [x*y, z].
        x,y dims correspond to pixel coordinates for each image
        z dim corresponds to hyperspectral image wavelength.
    '''
    imgs = [i for i in os.listdir(hyp_img_dir) if i.endswith('.png')]
    sorted_imgs = sorted(imgs, key=lambda x: int(x.split('_')[0]))
    all_arrs = []
    for i in sorted_imgs[2:]:
        img = np.array(Image.open('%s/%s' % (hyp_img_dir, i)).convert('L'))
        all_arrs.append(img)
    arrs = np.stack(all_arrs, axis=2)
    np.save('%s.npy' % hyp_img_dir, arrs)


def three_to_two(input_array_fn,
                 crop_dim=(1, 80, 320, 479),
                 output_format='csv'):
    '''
    convert 3d npy to 2d
    '''
    npy = np.load(input_array_fn)
    if crop_dim:
        left, up, right, down = crop_dim
        npy = npy[int(up):int(down),int(left):int(right), :]
    h, w, d = npy.shape
    npy_2d = npy.reshape(h*w, d)
    opp = input_array_fn.rsplit('.', 1)[0]
    if output_format == 'csv':
        out_fn = f"{opp}.2d.csv"
        np.savetxt(out_fn, npy_2d, delimiter=",")
    else:
        out_fn = f"{opp}.2d.npy"
        np.save(out_fn, npy_2d.astype(np.float64))
    print('Done!')


def inference(pretrained_model, input_np_array_fn):
    """
    using trained model to make predictions on 2d or 3d numpy array files. The
    output is a numpy array object which has the same number of columns as
    the training data.
    """
    model = load_model(pretrained_model)
    test_npy = np.load(input_np_array_fn)
    npy_shape = test_npy.shape
    np_dim = len(npy_shape)
    test_npy_2d = test_npy.reshape(npy_shape[0]*npy_shape[1], npy_shape[2]) \
        if np_dim == 3 else test_npy
    pred_prob = model.predict(test_npy_2d)
    predictions = pred_prob.argmax(axis=1)
    opp = input_np_array_fn.split('/')[-1].split('.npy')[0]+'.prd'
    if np_dim == 3:
        predictions = predictions.reshape(npy_shape[0], npy_shape[1])
        df = pd.DataFrame(predictions)
        # 0: background; 1: leaf; 2: stem; 3: panicle
        df1 = df.replace(0, 255).replace(1, 127).replace(2, 253)\
            .replace(3, 190)
        df2 = df.replace(0, 255).replace(1, 201).replace(2, 192)\
            .replace(3, 174)
        df3 = df.replace(0, 255).replace(1, 127).replace(2, 134)\
            .replace(3, 212)
        arr = np.stack([df1.values, df2.values, df3.values], axis=2)
        sm.imsave(f'{opp}.png', arr)
    elif np_dim == 2:
        np.savetxt(f'{opp}.csv', predictions)
    else:
        raise Exception('accept either 2 or 3 dim numpy array!')


def leaft_counting_stat(infer_csv, opp, agreement_threshold=0.5,
                        lc_range='5-10'):
    """
    Calculate leaf counting metrics including CountDiff, AbsCountDiff, MSE,
    Agreement, r2, p_value. Make scatter, bar plots

    Arguments
    ---------
    infer_csv:
        csv file including ground truth and predicted results
    opp:
        prefix of output files
    agreement_threshold:
        value above threshold will be considered as agreement
    lc_range:
        range of leave counts counted for stat
    """
    df_infer = pd.read_csv(infer_csv)
    slope, intercept, r_value, p_value, std_err = linregress(
        df_infer['groundtruth'], df_infer['prediction'])
    df_infer['prediction_adjusted'] = (df_infer['prediction']-intercept)/slope

    df_infer['diff'] = df_infer['groundtruth'] -\
        df_infer['prediction_adjusted']
    df_infer['abs_diff'] = np.abs(df_infer['diff'])
    mi, ma = df_infer['diff'].min(), df_infer['diff'].max()
    mi_int = np.ceil(mi) if mi > 0 else np.floor(mi)
    ma_int = np.ceil(ma) if ma > 0 else np.floor(ma)
    bins = np.arange(mi_int, ma_int+1)
    cats = pd.cut(df_infer['diff'], bins)
    ax1 = pd.value_counts(cats).sort_index().plot.bar(color='blue')
    plt.xticks(rotation=40)
    ax1.set_xlabel('Relative Count Differece')
    ax1.set_ylabel('Frequency')
    plt.tight_layout()
    plt.savefig(f'{opp}_diff_bar.png')
    print('diff bar/histogram plot done!')

    n = (df_infer['abs_diff'] <= float(agreement_threshold)).sum()
    agreement_ratio = n/df_infer.shape[0]
    rmse = np.sqrt(mean_squared_error(df_infer['groundtruth'],
                                      df_infer['prediction_adjusted']))
    bt = int(lc_range.split('-')[0])
    tp = int(lc_range.split('-')[1])
    x = np.array([bt - 0.5, tp + 0.5])
    y = slope*x+intercept
    mean, std = df_infer['diff'].mean(), df_infer['diff'].std()
    abs_mean, abs_std = df_infer['abs_diff'].mean(), df_infer['abs_diff'].std()
    txt = 'CountDiff: %.2f(%.2f)\n' % (mean, std)
    txt += 'AbsCountDiff: %.2f(%.2f)\n' % (abs_mean, abs_std)
    txt += 'r2: %.2f\n' % r_value**2
    txt += 'p value: %s\n' % p_value
    txt += 'RMSE: %.2f\n' % rmse
    txt += 'Agreement: %.2f\n' % agreement_ratio
    txt += 'y = %.2fx + %.2f' % (slope, intercept)
    with open('%s.statics' % opp, 'w') as f1:
        f1.write(txt)
    print('statistics done!')

    ax2 = df_infer.plot.scatter(x='groundtruth', y='prediction', alpha=0.5,
                                figsize=(7, 7), edgecolor='k')        
    ax2.plot(x, y, color='red', linewidth=2)
    ax2.text(x=bt, y=tp-2, s=txt, fontsize=12, color='red')
    plt.savefig('%s_scatter.png' % opp)
    print('scatter plot done!')
