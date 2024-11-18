from sklearn.metrics import roc_auc_score, f1_score, accuracy_score, recall_score, confusion_matrix

def calculate_auc_score(ground_truth, prediction):
    """
    Calculate the Area Under the Receiver Operating Characteristic Curve (AUC) score.

    *AUC is a common metric used to evaluate the performance of a binary classification model.
    It measures the area under the ROC curve, which plots the true positive rate against
    the false positive rate at various classification thresholds.

    Parameters
    ----------
    ground_truth : array-like of shape (n_samples,)
        True binary labels.
    prediction : array-like of shape (n_samples,)
        Predicted probabilities or decision values for the positive class.

    Returns
    -------
    float
        The AUC score, a value between 0 and 1, with higher values indicating better
        classification performance.
    """
    auc_score = roc_auc_score(ground_truth, prediction)

    return auc_score

def calculate_f1_score(ground_truth, prediction):
    f1_score = f1_score(ground_truth, prediction)

    return f1_score

def calculate_accuracy_score(ground_truth, prediction):
    accuracy_score = accuracy_score(ground_truth, prediction)

    return accuracy_score

def calculate_sensitivity_score(ground_truth, prediction):
    sensitivity_score = recall_score(ground_truth, prediction)

    return sensitivity_score

def calculate_specificity_score(ground_truth, prediction):

    tn, fp, fn, tp = confusion_matrix(ground_truth, prediction).ravel()

    # Calculate specificity
    specificity_score = tn / (tn + fp)

    return specificity_score
