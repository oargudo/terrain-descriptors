#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "core.h"
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QFile>
#include <QTextStream>
#include <QPainter>
#include <QButtonGroup>
#include <QStandardItemModel>
#include <iostream>
#include <fstream>
#include <sstream>


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    widget = new TerrainWidget();
	QGridLayout* GLlayout = new QGridLayout;
    GLlayout->addWidget(widget, 0, 0);
	GLlayout->setContentsMargins(0, 0, 0, 0);
	ui->glWidget->setLayout(GLlayout);

    widget->setCamera(Camera(Vector3(-10.0, -10.0, 10.0), Vector3(0.0, 0.0, 0.0)));

    createActions();
    createPresets();
    updatePresetInfo();

    // group radio buttons in DEM toolbox to be mutually exclusive across pages
    baseTexButtonGroup = new QButtonGroup(this);
    QList<QRadioButton*> allRadioButtons = ui->demToolbox->findChildren<QRadioButton*>();
    for (QRadioButton* rb : allRadioButtons) {
        baseTexButtonGroup->addButton(rb);
    }
}

MainWindow::~MainWindow()
{
    delete baseTexButtonGroup;
    delete ui;
}

void MainWindow::createActions()
{
    connect(widget, SIGNAL(glInitialized()), this, SLOT(loadPreset())); // first call
    connect(widget, SIGNAL(rayQueried(const Ray&)), this, SLOT(queryRay(const Ray&)));

    connect(ui->actionLoadPNG, SIGNAL(triggered()), this, SLOT(loadPNG()));
    connect(ui->actionLoadASC, SIGNAL(triggered()), this, SLOT(loadASC()));
    connect(ui->btn_loadPreset, SIGNAL(clicked()), this, SLOT(loadPreset()));
    connect(ui->cb_preset, SIGNAL(currentIndexChanged(int)), this, SLOT(updatePresetInfo()));
    connect(ui->btn_updateDEM, SIGNAL(clicked()), this, SLOT(updateDEM()));

    connect(ui->actionSaveTexture, SIGNAL(triggered()), this, SLOT(saveLayersCombined()));
    connect(ui->actionSaveLayerMetric, SIGNAL(triggered()), this, SLOT(saveLayerMetric()));
    connect(ui->actionSaveLayerRivers, SIGNAL(triggered()), this, SLOT(saveLayerRivers()));
    connect(ui->actionSaveLayerRidges, SIGNAL(triggered()), this, SLOT(saveLayerRidges()));
    connect(ui->actionSaveMetric, SIGNAL(triggered()), this, SLOT(saveMetric()));

    connect(ui->btn_fillDepressions, SIGNAL(clicked()), this, SLOT(editFillDepressions()));
    connect(ui->btn_breaching, SIGNAL(clicked()), this, SLOT(editBreaching()));
    connect(ui->btn_gaussianSmooth, SIGNAL(clicked()), this, SLOT(editGaussianSmooth()));

    connect(ui->dem_shading, SIGNAL(clicked()), this, SLOT(updateTexture()));
    connect(ui->dem_shadows, SIGNAL(clicked()), this, SLOT(updateLightmap()));
    connect(ui->light_azimuth,  SIGNAL(valueChanged(double)), this, SLOT(updateLightmap()));
    connect(ui->light_altitude, SIGNAL(valueChanged(double)), this, SLOT(updateLightmap()));
    connect(ui->fixed_palette, SIGNAL(clicked()), this, SLOT(modifyPalette()));
    connect(ui->cb_palette, SIGNAL(currentIndexChanged(int)), this, SLOT(modifyPalette()));
    connect(ui->sb_paletteMin, SIGNAL(valueChanged(double)), this, SLOT(modifyPalette()));
    connect(ui->sb_paletteMax, SIGNAL(valueChanged(double)), this, SLOT(modifyPalette()));
    connect(ui->btn_resetPaletteRange, SIGNAL(clicked()), this, SLOT(resetPaletteRange()));

    connect(ui->dem_uniform, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_elevation, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_waterlevel, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_waterlevel_m, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_waterlevel_relief, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_normals, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_slopeGradient, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_slopeAvg, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_aspect, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_aspect_east, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_aspect_north, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_laplacian, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_fractLaplacian, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_asymmetry, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_asymmetryDir, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_asymmetryDirTol, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_asymmetryW, SIGNAL(valueChanged(int)), this, SLOT(redrawBaseTexture()));

    connect(ui->dem_curv_fromQuadrics, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_curv_fromQuadrics_w, SIGNAL(valueChanged(int)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_curvMin, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_curvMax, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_curvMean, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_curvGaussian, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_curvContour, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_curvTangent, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_curvProfile, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));

    connect(ui->dem_localRelief_w, SIGNAL(valueChanged(int)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_localRelief, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_localVariance, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_areaRatio, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_tpi, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_ruggedness, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_surfaceRoughness, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_surfaceRoughnessStd, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));

    connect(ui->btn_pickViewshedLocation, SIGNAL(clicked()), this, SLOT(updateViewshedLocation()));
    connect(ui->dem_viewshedLoc, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_viewshedTotalIn, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_viewshedTotalOut, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_posOpenness, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_negOpenness, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_accessibility, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_skyview, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_skyviewApprox, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_sunlight, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_dahi, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_dahiMax, SIGNAL(valueChanged(int)), this, SLOT(redrawBaseTexture()));

    connect(ui->dem_showSingleLandform, SIGNAL(toggled(bool)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_showSingleLandform_id, SIGNAL(valueChanged(int)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_dikauWoodClass, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_fisher2004, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_fisher2004entropy, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_tpiClass, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_geomorphons, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_tophat, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));

    connect(ui->dem_streamArea, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamAreaLog, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamArea_exp, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamArea_flowAlgorithm, SIGNAL(currentIndexChanged(int)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamArea_mfdExp, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamArea_rho8iters, SIGNAL(valueChanged(int)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamArea_uniqueSource, SIGNAL(toggled(bool)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamArea_source_x, SIGNAL(valueChanged(int)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamArea_source_y, SIGNAL(valueChanged(int)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_wetness, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamPower, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamPower_m, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamPower_n, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamPower_threshold, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_streamPower_t, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseTexture()));    
    connect(ui->dem_depthToWater, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_hand, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_riverDistEuclidean, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_riverDistFlow, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_branchLength, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_branchLength_exp, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseTexture()));

    connect(ui->dem_peakedness, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_ors, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_jut, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_rut, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_jut_curvature, SIGNAL(toggled(bool)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_zreduced_converging, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_zreduced_diverging, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_zreduced_x, SIGNAL(valueChanged(int)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_zreduced_y, SIGNAL(valueChanged(int)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_deng2008_peakProto, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_deng2008_peakSimil, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_deng2008_peakPercent, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_nni, SIGNAL(clicked()), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_nni_radius, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseTexture()));
    connect(ui->dem_nni_relative, SIGNAL(toggled(bool)), this, SLOT(redrawBaseTexture()));

    connect(ui->showRivers, SIGNAL(clicked()), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_cit_s, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseRiversLayer()));
    connect(ui->rivers_cit_t, SIGNAL(valueChanged(double)), this, SLOT(redrawBaseRiversLayer()));
    connect(ui->rivers_scale, SIGNAL(valueChanged(double)), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_color, SIGNAL(currentIndexChanged(int)), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_width, SIGNAL(currentIndexChanged(int)), this, SLOT(redrawRiversLayer()));
    connect(ui->showSegmentSmooth, SIGNAL(clicked()), this, SLOT(redrawRiversLayer()));
    connect(ui->showRiversSimple, SIGNAL(clicked()), this, SLOT(redrawRiversLayer()));
    connect(ui->showRiversSources, SIGNAL(clicked()), this, SLOT(redrawRiversLayer()));
    connect(ui->showRiversJunctions, SIGNAL(clicked()), this, SLOT(redrawRiversLayer()));

    connect(ui->showBasins, SIGNAL(clicked()), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_drainageBasin_perimeter, SIGNAL(clicked()), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_drainageBasin_color, SIGNAL(currentIndexChanged(int)), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_drainageBasinAll, SIGNAL(clicked()), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_drainageBasinStrahler, SIGNAL(clicked()), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_drainageBasinStrahler_order, SIGNAL(valueChanged(int)), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_drainageBasinCoords, SIGNAL(clicked()), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_basin_x, SIGNAL(valueChanged(int)), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_basin_y, SIGNAL(valueChanged(int)), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_basin_useRiver, SIGNAL(clicked()), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_stopSea, SIGNAL(clicked()), this, SLOT(redrawRiversLayer()));
    connect(ui->rivers_stopSea_elev, SIGNAL(valueChanged(double)), this, SLOT(redrawRiversLayer()));

    connect(ui->showRidges, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_drainage, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_laplacian, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_curvature, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_drainageThreshold, SIGNAL(valueChanged(double)), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_laplacianThreshold, SIGNAL(valueChanged(double)), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_curvatureThreshold, SIGNAL(valueChanged(double)), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_PPA, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_PPA_profileLength, SIGNAL(valueChanged(int)), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_showPPAdebug, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_divtree, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_divtree_promCutoff, SIGNAL(valueChanged(int)), this, SLOT(redrawBaseRidgesLayer()));
    connect(ui->ridges_divtree_detailed, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_level, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_level_type, SIGNAL(currentIndexChanged(int)), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_level_scale, SIGNAL(valueChanged(double)), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_showPeaks, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_showSaddles, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_ignoreSea, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_ignoreSea_elev, SIGNAL(valueChanged(double)), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_showProm, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_displayPeakId, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_showIsol, SIGNAL(clicked()), this, SLOT(redrawRidgesLayer()));
    connect(ui->ridges_showIsol_peak, SIGNAL(valueChanged(int)), this, SLOT(redrawRidgesLayer()));
}

void MainWindow::updateHeightfield(bool resetCam)
{    
    if (terrainAnalysis) delete terrainAnalysis;
    terrainAnalysis = new TerrainAnalysis(hf);

    widget->setHeightfield(hf, resetCam);

    computeLightmap(ui->dem_shadows->isChecked());
    computeBaseTexture();
    computeRiversLayer();
    computeRidgesLayer();
    updateTexture();

    displayHeightfieldDimensions();

    cursorCell = Index2(hf.getSizeX()/2, hf.getSizeY()/2);
    ui->viewshed_location_x->blockSignals(true);
    ui->viewshed_location_y->blockSignals(true);
    ui->viewshed_location_x->setMaximum(hf.getSizeX());
    ui->viewshed_location_y->setMaximum(hf.getSizeY());
    ui->viewshed_location_x->setValue(cursorCell.x());
    ui->viewshed_location_y->setValue(cursorCell.y());
    ui->viewshed_location_x->blockSignals(false);
    ui->viewshed_location_y->blockSignals(false);

    ui->dem_streamArea_source_x->blockSignals(true);
    ui->dem_streamArea_source_y->blockSignals(true);
    ui->dem_streamArea_source_x->setMaximum(hf.getSizeX());
    ui->dem_streamArea_source_y->setMaximum(hf.getSizeY());
    ui->dem_streamArea_source_x->setValue(cursorCell.x());
    ui->dem_streamArea_source_y->setValue(cursorCell.y());
    ui->dem_streamArea_source_x->blockSignals(false);
    ui->dem_streamArea_source_y->blockSignals(false);

    ui->dem_zreduced_x->blockSignals(true);
    ui->dem_zreduced_y->blockSignals(true);
    ui->dem_zreduced_x->setMaximum(hf.getSizeX());
    ui->dem_zreduced_y->setMaximum(hf.getSizeY());
    ui->dem_zreduced_x->setValue(cursorCell.x());
    ui->dem_zreduced_y->setValue(cursorCell.y());
    ui->dem_zreduced_x->blockSignals(false);
    ui->dem_zreduced_y->blockSignals(false);

    ui->rivers_basin_x->blockSignals(true);
    ui->rivers_basin_y->blockSignals(true);
    ui->rivers_basin_x->setMaximum(hf.getSizeX());
    ui->rivers_basin_y->setMaximum(hf.getSizeY());
    ui->rivers_basin_x->setValue(cursorCell.x());
    ui->rivers_basin_y->setValue(cursorCell.y());
    ui->rivers_basin_x->blockSignals(false);
    ui->rivers_basin_y->blockSignals(false);
}

void MainWindow::updateLightmap()
{
    computeLightmap(ui->dem_shadows->isChecked());
    updateTexture();
}

void MainWindow::updateTexture()
{
    // update viewer
    widget->setTexture(composeTexture());
    widget->update();
}

void MainWindow::redrawBaseTexture()
{
    computeBaseTexture();
    updateTexture();
}

void MainWindow::redrawRiversLayer()
{
    computeRiversLayer();
    updateTexture();
}

void MainWindow::redrawRidgesLayer()
{
    computeRidgesLayer();
    updateTexture();
}

void MainWindow::redrawBaseRiversLayer()
{
    computeRiversLayer();
    if (ui->dem_depthToWater->isChecked() ||
        ui->dem_hand->isChecked() ||
        ui->dem_riverDistEuclidean->isChecked() ||
        ui->dem_riverDistFlow->isChecked())
    {
        computeBaseTexture();
    }
    updateTexture();
}

void MainWindow::redrawBaseRidgesLayer()
{
    computeRidgesLayer();
    if (ui->dem_nni->isChecked())
        computeBaseTexture();
    updateTexture();
}

void MainWindow::modifyPalette()
{
    computeBaseTexture(false);
    updateTexture();
}

void MainWindow::resetPaletteRange()
{
    setValueNoSignal(ui->sb_paletteMin, paletteMin / paletteScale);
    setValueNoSignal(ui->sb_paletteMax, paletteMax / paletteScale);
    computeBaseTexture(false);
    updateTexture();
}

void MainWindow::setStatusText(const QString& s)
{
    statusBar()->showMessage(s);
}


void MainWindow::resizeHeightfield()
{
    int newSizeX = ui->hf_gridX->value();
    int newSizeY = ui->hf_gridY->value();
    double minZ = ui->hf_elevMin->value();
    double maxZ = ui->hf_elevMax->value();
    double kmX = ui->hf_kmX->value()*1000;
    double kmY = ui->hf_kmY->value()*1000;

    HeightField aux = hf.setResolution(ui->hf_gridX->value(), ui->hf_gridY->value());
    double rmin, rmax;
    aux.getRange(rmin, rmax);
    std::vector<double> values(newSizeX * newSizeY);
    for (int i = 0; i < newSizeX * newSizeY; i++) {
        values[i] = (maxZ - minZ) * (aux[i] - rmin) / (rmax - rmin) + minZ;
    }

    hf = HeightField(Box2(kmX, kmY), newSizeX, newSizeY, values);
    updateHeightfield(true);
}

void MainWindow::editFillDepressions()
{
    hf.fillDepressions(0.1);
    updateHeightfield(false);
}

void MainWindow::editBreaching()
{
    hf.completeBreach();
    updateHeightfield(false);
}

void MainWindow::editGaussianSmooth()
{
    hf = HeightField(hf.gaussianBlur(ui->sb_gaussianBlur_radius->value()));
    updateHeightfield(false);
}

QImage MainWindow::composeTexture() const
{
    QImage background = baseTexture;
    if (ui->dem_shading->isChecked()) {
        background = shadeMap(background, shadeType, ui->dem_shadows->isChecked());
    }

    int upscale = std::max(scaleLayerRidges, scaleLayerRivers);
    QImage tex = background.scaled(upscale * background.width(), upscale * background.height(),
                    Qt::KeepAspectRatio, Qt::SmoothTransformation);
    QPainter painter(&tex);
    if (!riversLayer.isNull())
        painter.drawImage(QPoint(0, 0), riversLayer.scaled(tex.size(), Qt::IgnoreAspectRatio, Qt::SmoothTransformation));
    if (!ridgesLayer.isNull())
        painter.drawImage(QPoint(0, 0), ridgesLayer.scaled(tex.size(), Qt::IgnoreAspectRatio, Qt::SmoothTransformation));

    return tex;
}

void MainWindow::displayHeightfieldDimensions()
{
    ui->hf_gridX->setValue(hf.getSizeX());
    ui->hf_gridY->setValue(hf.getSizeY());

    Box3 box = hf.getBox();
    ui->hf_kmX->setValue(box.size()[0] / 1000.0);
    ui->hf_kmY->setValue(box.size()[1] / 1000.0);
    ui->hf_elevMin->setValue(box.getMin()[2]);
    ui->hf_elevMax->setValue(box.getMax()[2]);
    ui->label_cellSize->setText(QString::number(hf.getCellSize()[0], 'f', 1) + " x " +
            QString::number(hf.getCellSize()[1], 'f', 1) + " m");
}

void MainWindow::updatePaletteRangeSelectors()
{
    double vmin, vmax;
    currentMetric.getRange(vmin, vmax);
    vmin = std::min(vmin, paletteMin);
    vmax = std::max(vmax, paletteMax);

    double absLargest = std::max(std::abs(paletteMin), std::abs(paletteMax));
    if (absLargest < 0.1 || absLargest >= 100.0) {
        int pow10 = int(std::floor(std::log10(absLargest)));
        paletteScale = std::pow(10, pow10);
        ui->labelPaletteRange->setText(QStringLiteral("Range (e%1)").arg(pow10));
    }
    else {
        paletteScale = 1;
        ui->labelPaletteRange->setText("Range");
    }

    ui->sb_paletteMin->blockSignals(true);
    ui->sb_paletteMax->blockSignals(true);

    ui->sb_paletteMin->setMinimum(vmin / paletteScale);
    ui->sb_paletteMin->setMaximum(vmax / paletteScale);
    ui->sb_paletteMax->setMinimum(vmin / paletteScale);
    ui->sb_paletteMax->setMaximum(vmax / paletteScale);
    ui->sb_paletteMin->setValue(paletteMin / paletteScale);
    ui->sb_paletteMax->setValue(paletteMax / paletteScale);

    ui->sb_paletteMin->blockSignals(false);
    ui->sb_paletteMax->blockSignals(false);
}

void MainWindow::updateCursorInfoTexts()
{
    if (cursorCell.x() < 0 || cursorCell.y() < 0) return;

    double h = hf.at(cursorCell);
    ui->lineClickInfoCell->setText("(" + QString::number(cursorCell.x()) + ", " +
                                   QString::number(cursorCell.y()) + ")");
    ui->lineClickInfoElev->setText(QString::number(h, 'f', 2) + " m");
    ui->lineClickInfoMetric->setText(QString::number(currentMetric.at(cursorCell)));
}

void MainWindow::setValueNoSignal(QSpinBox *s, int val)
{
    s->blockSignals(true);
    s->setValue(val);
    s->blockSignals(false);
}

void MainWindow::setValueNoSignal(QDoubleSpinBox *s, double val)
{
    s->blockSignals(true);
    s->setValue(val);
    s->blockSignals(false);
}

void MainWindow::updateDEM()
{
    if (hf.getNumElements() == 0) return;
    resizeHeightfield();
}

void MainWindow::loadHeightfield(const QString& hfPath, double hfWidth, double hfHeight, double hfMinZ, double hfMaxZ, double scale)
{
    std::cout << "Loading " << hfPath.toStdString() << std::endl;
    QImage imghf = QImage(hfPath).mirrored(false, true);
    if (imghf.isNull()) {
        std::cout << "! Could not find: " << QFileInfo(hfPath).absoluteFilePath().toStdString() << std::endl;
        return;
    }

    imghf = imghf.scaledToWidth(int(scale * imghf.width() + 0.5), Qt::SmoothTransformation);

    hf = HeightField(Box2(hfWidth, hfHeight), imghf);
    hf.normalize();
    hf *= (hfMaxZ - hfMinZ);
    hf += (hfMinZ);

    updateHeightfield();
}

void MainWindow::loadPNG()
{
    QString filename = QFileDialog::getOpenFileName(nullptr,
                       "Choose a filename to load", "", "png file (*.png)");
    if (filename.isNull())
        return;

    QImage img(filename);
    hf = HeightField(Box2(ui->hf_kmX->value() * 1000.0, ui->hf_kmY->value() * 1000.0),
                     img.mirrored(false, true), ui->hf_elevMin->value(), ui->hf_elevMax->value());
    updateHeightfield();
}

void MainWindow::loadASC()
{
    QString filename = QFileDialog::getOpenFileName(nullptr,
                       "Choose a filename to load", "", "Ascii grid file (*.asc)");

    if (filename.isNull())
        return;

    std::string path = filename.toStdString();
    std::fstream infile(path, std::fstream::in);

    // read headers until we find a line starting with a float value
    int ncols = 0, nrows = 0;
    float minLat = 0, minLon = 0, cellSize = 0, nodataval = -9999;
    std::string line;
    bool centered = false;

    while (std::getline(infile, line)) {

        // get first "word" in the line
        std::istringstream liness(line);
        std::string s;
        liness >> s;

        // check if it is a float, in that case we finished headers
        std::istringstream iss(s);
        float f;
        if ((iss >> std::ws >> f) && !iss.fail()) {
            // go read the values
            break;
        }

        // otherwise, read value and save to appropiate variable
        for (unsigned int i = 0; i < s.size(); i++)
            s[i] = std::tolower(s[i]);

        if (s == "ncols")
            liness >> ncols;
        else if (s == "nrows")
            liness >> nrows;
        else if (s == "xllcorner")
            liness >> minLon;
        else if (s == "yllcorner")
            liness >> minLat;
        else if (s == "xllcenter") {
            liness >> minLon;
            centered = true;
        }
        else if (s == "yllcenter") {
            liness >> minLat;
            centered = true;
        }
        else if (s == "cellsize")
            liness >> cellSize;
        else if (s == "nodata_value")
            liness >> nodataval;
        else
            std::cerr << "Unrecognized header string " << s << std::endl;
    }
    if (centered) {
        minLat -= 0.5*cellSize;
        minLon -= 0.5*cellSize;
    }

    // now read the elevation values
    int num_samples = nrows * ncols;
    std::vector<double> samples(num_samples);
    int i = 0;
    std::istringstream iss(line);
    double val;
    while (iss >> val) {
        if (val > nodataval) samples[i] = val;
        else                 samples[i] = 0; // choose another elev?
        i++;
    }
    for (; i < num_samples; i++) {
        infile >> val;
        if (val > nodataval) samples[i] = val;
        else                 samples[i] = 0; // choose another elev?
    }
    infile.close();

    // flip DEM vertically
    hf = HeightField(HeightField(Box2((ncols - 1) * cellSize, (nrows - 1) * cellSize), ncols, nrows));
    for (int i = 0; i < ncols; i++) {
        for (int j = 0; j < nrows; j++) {
            hf(i, j) = samples[hf.cellId(i, nrows - 1 - j)];
        }
    }

    updateHeightfield();
}


void MainWindow::saveTexture(const QImage& img) const
{
    QString filename = QFileDialog::getSaveFileName(nullptr,
                                    "Choose a filename to save",
                                    "",
                                    "PNG file (*.png)");
    if (filename.isNull())
        return;

    img.mirrored(false, true).save(filename);
}

void MainWindow::saveMetric()
{
    QString filename = QFileDialog::getSaveFileName(nullptr,
                                    "Choose a filename to save",
                                    "",
                                    "TXT file (*.txt)");
    if (filename.isNull())
        return;

    QFile file(filename);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "Failed to open file for writing";
        return;
    }
    QTextStream out(&file);
    for (int i = 0; i < currentMetric.getSizeX(); i++) {
        for (int j = 0; j < currentMetric.getSizeY(); j++) {
            if (j > 0) out << " ";
            out << currentMetric.at(i, j);
        }
        out << "\n";
    }
    file.close();
}


void MainWindow::loadPreset()
{
    if (ui->cb_preset->count() == 0) return;

    std::string selected = ui->cb_preset->currentText().toStdString();
    if (presetTerrains.find(selected) == presetTerrains.end()) return;

    const PresetTerrain& preset = presetTerrains[selected];
    loadHeightfield(preset.imagePath, preset.terrainX, preset.terrainY,
                    preset.hmin, preset.hmax, ui->sb_presetFactor->value());
    ui->dem_sunlightLatitude->setValue(preset.avgLat);
}

void MainWindow::updatePresetInfo()
{
    const PresetTerrain& preset = presetTerrains[ui->cb_preset->currentText().toStdString()];
    ui->label_preset->setText(
                QString::number(preset.cellsX) + " x "
                + QString::number(preset.cellsY) + " cells,  "
                + QString::number(preset.terrainX/1000.0, 'f', 1) + " x "
                + QString::number(preset.terrainY/1000.0, 'f', 1) + " km"
                );
    ui->sb_presetFactor->setValue(preset.defaultScale);
}

void disableComboBoxItem(QComboBox * comboBox, int index)
{
    QStandardItemModel* model = qobject_cast<QStandardItemModel*>(comboBox->model());
    if (!model) return;

    QStandardItem* item = model->item(index);
    if (!item) return;
    item->setEnabled(false);

    QFont f = item->font();
    f.setBold(true);
    item->setFont(f);
}

void MainWindow::createPresets()
{
    QFile file("terrains/presets.txt");
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        std::cout << "Could not find preset terrains file" << std::endl;
        return;
    }

    ui->cb_preset->clear();
    presetTerrains.clear();

    QTextStream in(&file);
    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        bool isSelected = false;

        if (line.isEmpty()) {
            ui->cb_preset->addItem("");
            disableComboBoxItem(ui->cb_preset, ui->cb_preset->count() -1);
            continue;
        }
        if (line.startsWith("#")) {
            QString header = line.mid(1).trimmed();
            ui->cb_preset->addItem(header);
            disableComboBoxItem(ui->cb_preset, ui->cb_preset->count() - 1);
            continue;
        }
        if (line.startsWith("!")) {
            isSelected = true;
            line = line.mid(1);
        }

        QStringList parts = line.split(',', Qt::SkipEmptyParts);
        if (parts.size() != 10)
            continue;

        QString name = parts[0].trimmed().remove('\"');
        QString imagePath = parts[1].trimmed().remove('\"');
        int imgWidth = parts[2].toInt();
        int imgHeight = parts[3].toInt();
        int terrainWidth = parts[4].toInt();
        int terrainHeight = parts[5].toInt();
        double minElev = parts[6].toDouble();
        double maxElev = parts[7].toDouble();
        double latitude = parts[8].toDouble();
        double scale = parts[9].toDouble();

        ui->cb_preset->addItem(name);
        if (isSelected)
            ui->cb_preset->setCurrentIndex(ui->cb_preset->count()-1);
        presetTerrains[name.toStdString()] = PresetTerrain(
            dataPath + imagePath,
            imgWidth, imgHeight,
            terrainWidth, terrainHeight,
            minElev, maxElev,
            latitude, scale
        );
    }
}


void MainWindow::queryRay(const Ray& ray)
{
    double t;
    Vector3 p;
    if (hf.Intersect(ray, t, p, hf.getBox(), 1)) {
        Vector3 poff = p + Vector3(0.5*hf.getCellSize()[0], 0.5*hf.getCellSize()[1], 0);
        cursorCell = hf.cellCoords(Vector2(poff));

        widget->setCursorStyle(Norm(hf.getCellSize()), Vector3(0.0, 0.7, 0));
        widget->showCursor(p);
        updateCursorInfoTexts();

        if (ui->btn_pickStreamAreaSource->isChecked()) {
            setValueNoSignal(ui->dem_streamArea_source_x, cursorCell.x());
            setValueNoSignal(ui->dem_streamArea_source_y, cursorCell.y());
            if (ui->dem_streamArea_uniqueSource->isChecked() &&
                (ui->dem_streamArea->isChecked()  || ui->dem_streamAreaLog->isChecked() ||
                 ui->dem_streamPower->isChecked() || ui->dem_wetness->isChecked())) {
                redrawBaseTexture();
            }
        }
        if (ui->btn_pickZReduced->isChecked()) {
            setValueNoSignal(ui->dem_zreduced_x, cursorCell.x());
            setValueNoSignal(ui->dem_zreduced_y, cursorCell.y());
            if (ui->dem_zreduced_converging->isChecked() || ui->dem_zreduced_diverging->isChecked()) {
                redrawBaseTexture();
            }
        }
        if (ui->btn_pickBasinLocation->isChecked()) {
            setValueNoSignal(ui->rivers_basin_x, cursorCell.x());
            setValueNoSignal(ui->rivers_basin_y, cursorCell.y());
            if (ui->showBasins->isChecked() && ui->rivers_drainageBasinCoords->isChecked()) {
                redrawRiversLayer();
            }
        }
    }
    else {
        cursorCell = Index2(-1,-1);
        widget->hideCursor();
        ui->lineClickInfoCell->setText("");
        ui->lineClickInfoElev->setText("");
        ui->lineClickInfoMetric->setText("");
    }
}

void MainWindow::updateViewshedLocation()
{
    ui->viewshed_location_x->setValue(cursorCell.x());
    ui->viewshed_location_y->setValue(cursorCell.y());
    if (ui->dem_viewshedLoc->isChecked()) {
        redrawBaseTexture();
    }
}
