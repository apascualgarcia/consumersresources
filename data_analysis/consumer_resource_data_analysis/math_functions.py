def closest_element_in_list(el, liste):
    closest_index=0
    distance=abs(el-liste[0])
    for i in range(1, len(liste)):
        dist = abs(el-liste[i])
        if(dist<=distance):
            distance=dist
            closest_index=i
    return liste[closest_index]

def index_closest_element_in_list(el, liste):
    closest_index=0
    distance=abs(el-liste[0])
    for i in range(1, len(liste)):
        dist = abs(el-liste[i])
        if(dist<=distance):
            distance=dist
            closest_index=i
    return closest_index
