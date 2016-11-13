// Catherine Galkina, group 524, year 2016
// File point.h

class Point {
    float coord[2];
    int index;
public:
    Point() { coord[0] = coord[1] = index = 0; }
    Point(float x, float y, int idx) : index(idx)
        { coord[0] = x; coord[1] = y; }
    float GetX() const { return coord[0]; }
    float GetY() const { return coord[1]; }
    int GetIndex() const { return index; }
};
