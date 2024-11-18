from scipy.stats import pearsonr

def calculate_overlapped_elements(set1, set2):
    """
    Calculate the number of overlapped elements between two sets.

    Parameters
    ----------
    set1 : set
        The first set.
    set2 : set
        The second set.

    Returns
    -------
    int
        The number of elements that are common to both sets.
    """
    return len(set1.intersection(set2))

def calculate_jaccard_coefficient(set1, set2):
    """
    Calculate the Jaccard coefficient between two sets.

    The Jaccard coefficient is a measure of similarity between two sets. It is
    defined as the size of the intersection divided by the size of the union
    of the sets.

    Parameters
    ----------
    set1 : set
        The first set.
    set2 : set
        The second set.

    Returns
    -------
    float
        The Jaccard coefficient between the two sets. Returns 0 if the union
        of sets is empty to avoid division by zero.
    """
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    coefficient = intersection / union if union != 0 else 0  # Avoid division by zero

    return coefficient

def calculate_dissimilarity(set1, set2):
    """
    Calculate the dissimilarity between two sets based on Pearson correlation coefficient.

    The dissimilarity is computed as 1 minus the Pearson correlation coefficient
    between the two sets.

    Parameters
    ----------
    set1 : array-like
        The first set or array for Pearson correlation calculation.
    set2 : array-like
        The second set or array for Pearson correlation calculation.

    Returns
    -------
    float
        The dissimilarity between the two sets, computed as 1 minus the Pearson
        correlation coefficient.
    """
    corr_coefficient, _ = pearsonr(set1, set2)

    dissimilarity = 1 - corr_coefficient

    return dissimilarity
