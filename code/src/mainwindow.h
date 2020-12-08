#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include "layer.h"
#include <QMainWindow>
#include <QPainter>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    MainWindow(const Layer* upperLayer = nullptr, const Layer* lowerLayer = nullptr,
               const GlobalParameters* gp = nullptr,
               const std::vector<Family*>& vecFamily = std::vector<Family*>(),
               QWidget *parent = 0);
    ~MainWindow();
    void setLayers(const Layer* upperLayer, const Layer* lowerLayer,
                   const std::vector<Family*>& vecFamily, const GlobalParameters* gp) {
        _upperLayer = upperLayer;
        _lowerLayer = lowerLayer;
        _vecFamily = vecFamily;
        _gp = gp;
    }
protected:
    void paintEvent(QPaintEvent *);    // 繪製背景圖
    void wheelEvent(QWheelEvent *);    // 滾輪縮放事件
    void keyPressEvent(QKeyEvent *event);	//鍵盤事件
    void mousePressEvent(QMouseEvent *event);	//滑鼠點下事件
    void mouseMoveEvent(QMouseEvent *event);	//滑鼠移動事件
    void resizeEvent(QResizeEvent *event);

private:
    void findMinBBox();
    void resetViewPoint();
    Ui::MainWindow *ui;
    const Layer *_upperLayer, *_lowerLayer;
    const GlobalParameters* _gp;
    std::vector<Family*> _vecFamily;
    int _win_width;      // 視窗寬度
    int _win_height;     // 視窗高度
    double _llx = DBL_MAX, _lly = DBL_MAX, _urx = DBL_MIN, _ury = DBL_MIN;
    double _scaleX = 1.0, _scaleY = 1.0;
    double _viewx = 0.0, _viewy = 0.0;
    int _numGridX = 30, _numGridY = 30, _minUnit = 0;
    bool _showSolution = false, _showText = false;
    double _base;
    bool _isResident = true, _isUpperExist = true, _isUpperCand = true, _isLowerExist = true,
    _isLowerCand = true, _isGrid = true, _isUpperSol = true, _isLowerSol = true, _isUpperNonSol = true, _isLowerNonSol = true;
    double _upperRatio = 0, _lowerRatio = 0, _famRatio = 0, _gridRatio = 0;
    bool _isLL = false, _isSavePic = false;
};

#endif // MAINWINDOW_H
