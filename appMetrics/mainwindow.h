#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QButtonGroup>
#include <QSpinBox>
#include "terrainwidget.h"
#include "heightfield.h"
#include "terrainanalysis.h"


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

public slots:
    void updateHeightfield(bool resetCam = true);
    void setStatusText(const QString& s);

    // update texture slots
    void updateLightmap();
    void updateTexture();
    void redrawBaseTexture();
    void redrawRiversLayer();
    void redrawRidgesLayer();
    void redrawBaseRiversLayer();
    void redrawBaseRidgesLayer();
    void modifyPalette();
    void resetPaletteRange();

    // terrain load actions
    void loadPNG();
    void loadASC();
    void updateDEM();

    // presets
    void loadPreset();
    void updatePresetInfo();

    // texture export actions
    void saveLayersCombined() { saveTexture(composeTexture()); };
    void saveLayerMetric() { saveTexture(baseTexture); };
    void saveLayerRivers() { saveTexture(riversLayer); };
    void saveLayerRidges() { saveTexture(ridgesLayer); };
    void saveMetric();

    // terrain edit actions
    void editFillDepressions();
    void editBreaching();
    void editGaussianSmooth();

    // queries
    void queryRay(const Ray& ray);

    // others
    void updateViewshedLocation();

private:
    void createActions();
    void createPresets();

    void loadHeightfield(const QString& hfPath, double hfWidth, double hfHeight, double hfMinZ, double hfMaxZ, double scale);
    void resizeHeightfield();

    void computeLightmap(bool doShadows);
    void computeMetric();
    void computeBaseTexture(bool updateMetric = true);
    void computeRiversLayer();
    void computeRidgesLayer();

    QImage composeTexture() const;
    QImage shadeMap(const QImage& img, int shadeType = 0, bool shadows = false) const;
    ColorPalette getPalette(const ColorPalette& preferred) const;

    void saveTexture(const QImage& img) const;

    void displayHeightfieldDimensions();
    void updatePaletteRangeSelectors();
    void updateCursorInfoTexts();
    void setValueNoSignal(QSpinBox* s, int val);
    void setValueNoSignal(QDoubleSpinBox* s, double val);


private:
    Ui::MainWindow *ui;
    TerrainWidget *widget;
    QButtonGroup *baseTexButtonGroup;
    QString dataPath = "./";

    HeightField hf;
    ScalarField2 currentMetric;
    TerrainAnalysis* terrainAnalysis = nullptr;
    Index2 cursorCell = Index2(-1,-1);

    Vector3 lightDirection;
    ScalarField2 lightMap, shadowMap;

    ColorPalette paletteDefault;
    double paletteMin = 0, paletteMax = 1;
    double paletteScale = 1;

    QImage baseTexture, riversLayer, ridgesLayer;
    int shadeType;
    const int scaleLayerRivers = 8;
    const int scaleLayerRidges = 8;
    static LookupPalette paletteLandDW, paletteLandTPI, paletteLandGeomorphons;
    static QVector<QString> namesLandDW, namesLandTPI, namesLandGeomorphons;

    struct PresetTerrain {
        QString imagePath = "";
        int cellsX = 0;
        int cellsY = 0;
        double terrainX = 0;
        double terrainY = 0;
        double hmin = 0;
        double hmax = 0;
        double avgLat = 45;
        double defaultScale = 1;

        PresetTerrain() {};
        PresetTerrain(const QString& img, int cx, int cy, double tx, double ty, double hmin, double hmax, double lat, double s)
                : imagePath(img), cellsX(cx), cellsY(cy), terrainX(tx), terrainY(ty),
                  hmin(hmin), hmax(hmax), avgLat(lat), defaultScale(s) {}
    };
    std::map<std::string, PresetTerrain> presetTerrains;
};
#endif // MAINWINDOW_H
