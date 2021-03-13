from .epijinn import Methylase


class Dnd(Methylase):
    pass


Dnd_EcoB7A = Dnd("Dnd_EcoB7A", "GAAC", 0, 3)
Dnd_Sli1326 = Dnd("Dnd_Sli1326", "GGCC", 0, 3)
Dnd_VciFF75 = Dnd("Dnd_VciFF75", "CCA", 0, 2)

DND = [Dnd_EcoB7A, Dnd_Sli1326, Dnd_VciFF75]
