#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QPainter>
#include <QWheelEvent>
#include <QtMath>

MainWindow::MainWindow(const Layer *upperLayer, const Layer *lowerLayer, const GlobalParameters* gp,
                       const std::vector<Family*>& vecFamily,
                       QWidget *parent) :
    QMainWindow(parent), ui(new Ui::MainWindow), _upperLayer(upperLayer), _lowerLayer(lowerLayer),
    _gp(gp), _vecFamily(vecFamily)
{
    ui->setupUi(this);
    _win_width = 1680;
    _win_height = 900;
    this->setFixedSize(_win_width, _win_height);	//視窗大小
    this->adjustSize();
    // this->setStyleSheet("background-color: white;");
    setWindowTitle(tr("Bi-level Competative Facility Locator"));

    findMinBBox();
    resetViewPoint();
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

void MainWindow::resetViewPoint(){
    _viewx = width()/2;
    _viewy = height()/2;
    _scaleX = _scaleY = (double) qMin(width(), height())/qMax(_urx, _ury)/2.2;
    _upperRatio = 0;
    _lowerRatio = 0;
    _gridRatio = 0;
    _famRatio = 0;
    _base = qSqrt(qAbs(_urx/30));
}

void MainWindow::findMinBBox(){
    if (_upperLayer)
        for (auto& fac : _upperLayer->vecFacility()) {
            _llx = qMin(fac->x(), _llx);
            _lly = qMin(fac->y(), _lly);
            _urx = qMax(fac->x(), _urx);
            _ury = qMax(fac->y(), _ury);
        }
    if (_lowerLayer)
        for (auto& fac : _lowerLayer->vecFacility()) {
            _llx = qMin(fac->x(), _llx);
            _lly = qMin(fac->y(), _lly);
            _urx = qMax(fac->x(), _urx);
            _ury = qMax(fac->y(), _ury);
        }

    if (!_vecFamily.empty())
        for (auto& fam : _vecFamily) {
            _llx = qMin(fam->x(), _llx);
            _lly = qMin(fam->y(), _lly);
            _urx = qMax(fam->x(), _urx);
            _ury = qMax(fam->y(), _ury);
        }

    _urx = qCeil((double)qMax(qAbs(_llx),qAbs(_urx))/_numGridX)*(1+_numGridX);
    _ury = qCeil((double)qMax(qAbs(_lly),qAbs(_ury))/_numGridY)*(1+_numGridY);
    _urx = _ury = qMax(_urx, _ury);
    // if (_llx < -1)
        _llx = -_urx;
    // if (_lly < -1)
        _lly = -_ury;
    _minUnit = qMin(_ury/_numGridY, _urx/_numGridX);
    _base = qSqrt(qAbs(_urx/30));
}

void MainWindow::resizeEvent(QResizeEvent *event)
{

    //scaled window size
    _win_width = double(event->size().width())/360.0;
    _win_height = double(event->size().height())/180.0;
}

void MainWindow::wheelEvent(QWheelEvent *event)
{
    // std::cout << "w _scaleX = " << _scaleX << std::endl;
    if (event->orientation() == Qt::Vertical){
        if (event->delta() > 0){
            _scaleX *= 1.01;
            _scaleY *= 1.01;
        } else {
            _scaleX /= 1.01;
            _scaleY /= 1.01;
        }
    }
    update();
}

void MainWindow::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        _viewx = event->x();
        _viewy = event->y();
        update();
    }
}

void MainWindow::mouseMoveEvent(QMouseEvent *event)
{
    _viewx = event->x();
    _viewy = event->y();
    update();
}

void MainWindow::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_F) {
        resetViewPoint();
    }
    if (event->key() == Qt::Key_T) {
        _showText ^= 1;
    }
    if (event->key() == Qt::Key_S) {
        _showSolution ^= 1;
    }
    if (event->key() == Qt::Key_Z) {
        _base += 0.01;
    }
    if (event->key() == Qt::Key_X) {
        if (_base > 0.01)
            _base -= 0.01;
    }
    if (event->key() == Qt::Key_1) {
        _isResident ^= 1;
    }
    if (event->key() == Qt::Key_2) {
        _isUpperCand ^= 1;
    }
    if (event->key() == Qt::Key_3) {
        _isUpperExist ^= 1;
    }
    if (event->key() == Qt::Key_4) {
        _isLowerCand ^= 1;
    }
    if (event->key() == Qt::Key_5) {
        _isLowerExist ^= 1;
    }
    if (event->key() == Qt::Key_6) {
        _isGrid ^= 1;
    }
    if (event->key() == Qt::Key_Q) {
        _isUpperSol ^= 1;
    }
    if (event->key() == Qt::Key_W) {
        _isLowerSol ^= 1;
    }
    if (event->key() == Qt::Key_E) {
        _isUpperNonSol ^= 1;
    }
    if (event->key() == Qt::Key_R) {
        _isLowerNonSol ^= 1;
    }
    if (event->key() == Qt::Key_7) {
        _famRatio += 0.1;
    }
    if (event->key() == Qt::Key_8) {
        _upperRatio += 0.1;
    }
    if (event->key() == Qt::Key_9) {
        _lowerRatio += 0.1;
    }
    if (event->key() == Qt::Key_0) {
        _gridRatio += 0.1;
    }
    if (event->key() == Qt::Key_L) {
        _viewx = width()*0.01;
        _viewy = height()*0.95;
    }
    if (event->key() == Qt::Key_P) {
        _isSavePic = true;
    }
    update();
}

// 繪圖事件, 用来產生背景
// 會自動一直 call 他
void MainWindow::paintEvent(QPaintEvent*)
{
    // setAttribute(Qt::WA_OpaquePaintEvent);
    QPainter painter(this);
    // this 代表 Main Window
    // 代表這個 painter 是要畫在 Main Window 上
    // 若要畫在一張圖上
    // 也可以自己創 QImage image; QPainter painter(&image);
    // QImage image;
    // QPainter painter(&image);
    // 畫完後 image.save(路徑, 圖檔類型);
    // clear old draw
    painter.eraseRect(0, 0, width(), height());
    painter.fillRect(0, 0, width(), height(), Qt::white);

    // put (0, 0) at center
    painter.translate(_viewx, _viewy);
    // painter.setTransform(QTransform::fromScale(-1, 1));
    // right = positive x, up = positive y

    // std::cout  << "base size = " << _base << std::endl;
    // scaling (zoom in/out)
    painter.scale(_scaleX, _scaleY);
    QFont font = painter.font();
    font.setPixelSize(qCeil(_base));
    painter.setFont(font);
    // painter->setStyleSheet("background-color: white;");
    // painter.setBackgroundMode(Qt::TransparentMode);
    // painter.setBackground(Qt::black);
    painter.setRenderHint(QPainter::Antialiasing, true);


    QTransform t = painter.transform();
    t.rotate(180, Qt::XAxis);
    painter.setTransform(t);


    painter.setPen(QPen(Qt::lightGray, _base/10/0.5*(1+_gridRatio), Qt::SolidLine));
    // draw grids
    if (_isGrid)
    for (int i = 1, j = _minUnit; j < _ury; ++i, j += _minUnit) {
        painter.drawLine(_llx, j, _urx, j);
        painter.drawLine(_llx, -j, _urx, -j);
    }
    if (_isGrid)
    for (int i = 1, j = _minUnit; j < _urx; ++i, j += _minUnit) {
        painter.drawLine(j, _lly, j, _ury);
        painter.drawLine(-j, _lly, -j, _ury);
    }
    // draw minimun bounding box
    // 筆是畫方塊外框，筆刷是畫方塊裡面
    painter.setPen(QPen(Qt::darkGray, _base/5*(1+_gridRatio), Qt::SolidLine)); // 設定筆：黑色、1.0px、實線
    // painter.setBrush(QBrush(Qt::gray, Qt::NoBrush)); // 設定筆刷：gray、square
    if (_isGrid){
    painter.drawRect(_llx, _lly, _urx-_llx, _ury-_lly);
    // draw x-axis, y-axis
    painter.drawLine(_llx, 0, _urx, 0);
    painter.drawLine(0, _lly, 0, _ury);
    }

    QPen pen;
    pen.setCapStyle(Qt::RoundCap);
    // double base = qSqrt(qAbs(_urx/20));

    if (_upperLayer) {
        pen.setColor(Qt::red);
        painter.setPen(pen);
        pen.setWidthF(_base);
        int i = 0;
        const PartSolution& par = _upperLayer->best_sol();
        pen.setWidthF(_base*(1+_upperRatio));
        if (_isUpperCand)
        for (auto& fac : _upperLayer->vecFacility()) {
            if (par.test(i) && _showSolution) {
                pen.setColor(Qt::magenta);
                painter.setPen(pen);
            } else {
                pen.setColor(Qt::red);
                painter.setPen(pen);
            }
            if (par.test(i) && !_isUpperSol) { i++; continue; }
            if (!par.test(i) && !_isUpperNonSol) { i++; continue; }
            QPointF pt (fac->x(), fac->y());
            painter.drawPoint(pt);
            i++;
        }

        pen.setColor(Qt::darkRed);
        painter.setPen(pen);
        if (_isUpperExist)
        for (auto& fac : _upperLayer->vecExistedFacility()) {
            QPointF pt (fac->x(), fac->y());
            painter.drawPoint(pt);
        }
    }
    if (_lowerLayer) {
        pen.setColor(Qt::cyan);
        painter.setPen(pen);
        pen.setWidthF(_base);
        int i = 0;
        const PartSolution& par = _lowerLayer->best_sol();
        pen.setWidthF(_base*(1+_lowerRatio));
        if (_isLowerCand)
        for (auto& fac : _lowerLayer->vecFacility()) {
            if (par.test(i) && _showSolution) {
                pen.setColor(Qt::blue);
                painter.setPen(pen);
            } else {
                pen.setColor(Qt::cyan);
                painter.setPen(pen);
            }
            if (par.test(i) && !_isLowerSol) { i++; continue; }
            if (!par.test(i) && !_isLowerNonSol) { i++; continue; }
            QPointF pt (fac->x(), fac->y());
            painter.drawPoint(pt);
            i++;
        }

        pen.setColor(Qt::darkCyan);
        painter.setPen(pen);
        if (_isLowerExist)
        for (auto& fac : _lowerLayer->vecExistedFacility()) {
            QPointF pt (fac->x(), fac->y());
            painter.drawPoint(pt);
        }


    }
    // std::cout << "maxbp = " << _gp->maxbp() << std::endl;
    // std::cout << "_vecFamily.empty() = " << _vecFamily.empty() << std::endl;
    if (!_vecFamily.empty() && _isResident) {
        pen.setColor(Qt::black);
        for (auto& fam : _vecFamily) {
            // std::cout << "fam bp: " << fam->bp() << ", maxbp = " << _gp->maxbp() << std::endl;
            // std::cout << "fam ratio: " << qAbs(fam->bp()/(double)_gp->maxbp()) << std::endl;
            pen.setWidthF((qAbs(fam->bp()/(double)_gp->maxbp())*_base/2+_base/2)*(1+_famRatio));
            painter.setPen(pen);
            QPointF pt (fam->x(), fam->y());
            painter.drawPoint(pt);
        }
    }

    if (_showText) {
        painter.setPen(QPen(Qt::black, _base/10.0, Qt::SolidLine));
        painter.drawText(0, _minUnit, "(0," + QString::number(_minUnit) + ")");
        painter.drawText(0, _ury, "(0," + QString::number(_ury) + ")");
        painter.drawText(_minUnit, 0, "(" + QString::number(_minUnit) + ",0)");
        painter.drawText(_urx, 0, "(" + QString::number(_urx) + ",0)");
    }

    painter.end();
    if (_isSavePic) {
        _isSavePic = false;
        // image.save((_gp->currentDateTime() + ".jpg").c_str(), "JPG");
    }
    /*
    // 畫字 (id)
    painter.drawText(rec[i].x+0.38*rec[i].w,rec[i].y+0.62*rec[i].h,
                     QString::number(rec[i].id));
    painter.setPen(QPen(Qt::red, 2, Qt::DashLine)); // 設定筆：紅色、2px、虛線
    painter.setBrush(Qt::NoBrush); // 設定筆刷：無筆刷
    // 畫邊界
    painter.drawRect(0,0,boundarySize,boundarySize);
    painter.end();
    */
}

MainWindow::~MainWindow()
{
    delete ui;
}
