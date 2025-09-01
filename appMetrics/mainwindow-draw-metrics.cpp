#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QPainter>
#include <iostream>
#include <queue>

#define ColorRGB(r, g, b) Vector3((r) / 255.0, (g) / 255.0, (b) / 255.0)

LookupPalette MainWindow::paletteLandDW({
    ColorRGB(0, 0, 0),
    ColorRGB(230, 151, 1),      // 1. nose
    ColorRGB(128, 164, 32),     // 2. shoulder slope
    ColorRGB(242, 91, 115),     // 3. hollow shoulder
    ColorRGB(255, 210, 128),    // 4. spur
    ColorRGB(161, 255, 143),    // 5. planar slope
    ColorRGB(255, 87, 255),     // 6. hollow
    ColorRGB(255, 255, 1),      // 7. spur foot
    ColorRGB(57, 168, 1),       // 8. slope foot
    ColorRGB(169, 1, 230),      // 9. hollow foot
    ColorRGB(255, 16, 1),       // 10. peak
    ColorRGB(255, 85, 1),       // 11. ridge
    ColorRGB(192, 192, 192),    // 12. plain
    ColorRGB(159, 252, 255),    // 13. saddle
    ColorRGB(1, 197, 255),      // 14. channel
    ColorRGB(1, 64, 255)        // 15. pit
});
QVector<QString> MainWindow::namesLandDW({
     "unclassified", "nose", "shoulder slope", "hollow shoulder",
     "spur", "planar slope", "hollow", "spur foot",
     "slope foot", "hollow foot", "peak", "ridge",
     "plain", "saddle", "channel", "pit"
});

LookupPalette MainWindow::paletteLandTPI({  // color scheme defined in [Weiss 2001]
    ColorRGB(192, 192, 192),
    ColorRGB(169, 1, 230),    // canyons, deeply incised streams
    ColorRGB(255, 87, 255),   // midslope drainages, shallow valleys
    ColorRGB(1, 112, 255),    // upland drainages, headwaters
    ColorRGB(1, 197, 255),    // U-shape valleys
    ColorRGB(57, 168, 1),     // plains
    ColorRGB(161, 255, 143),  // open slopes
    ColorRGB(255, 255, 1),    // upper slopes, mesas
    ColorRGB(255, 210, 128),  // local ridges/hills in valleys
    ColorRGB(230, 151, 1),    // midslope ridges, small hills in plains
    ColorRGB(255, 85, 1)      // peak
});
QVector<QString> MainWindow::namesLandTPI({
     "unclassified", "canyons, deeply incised streams", "midslope drainages, shallow valleys", "upland drainages, headwaters",
     "U-shape valleys", "plains", "open slopes", "upper slopes, mesas",
     "local ridges/hills in valleys", "midslope ridges, small hills in plains", "peaks, ridges"
});

LookupPalette MainWindow::paletteLandGeomorphons({
    ColorRGB(128, 128, 128), // 0. flat
    ColorRGB(255, 85, 1),    // 1. peak
    ColorRGB(230, 151, 1),   // 2. ridge
    ColorRGB(255, 210, 128), // 3. shoulder
    ColorRGB(161, 255, 143), // 4. hollow
    ColorRGB(57, 168, 1),    // 5. slope
    ColorRGB(255, 255, 1),   // 6. spur
    ColorRGB(1, 197, 255),   // 7. footslope
    ColorRGB(255, 87, 255),  // 8. valley
    ColorRGB(169, 1, 230),   // 9. pit
});
QVector<QString> MainWindow::namesLandGeomorphons({
     "flat", "peak", "ridge", "shoulder",
     "hollow", "slope", "spur", "footslope",
     "valley", "pit"
});


void MainWindow::computeLightmap(bool doShadows)
{
    double aAzimuth  = Math::DegreeToRadian(ui->light_azimuth->value());
    double aAltitude = Math::DegreeToRadian(ui->light_altitude->value());

    // note that azimuth is given CW from N
    lightDirection = Vector3(std::cos(aAltitude)*std::sin(aAzimuth),
                             std::cos(aAltitude)*std::cos(aAzimuth),
                             std::sin(aAltitude));

    lightMap  = hf.DirectLight(lightDirection);
    if (doShadows) {
        shadowMap = hf.SelfShadow(lightDirection);
    }
}

QImage MainWindow::shadeMap(const QImage& img, int shadeType, bool shadows) const
{
    QImage shading = img.copy();
    for (int i = 0; i < shading.width(); i++) {
        for (int j = 0; j < shading.height(); j++) {
            QColor qc = img.pixelColor(i, j);
            Vector3 cz = Vector3(qc.redF(), qc.greenF(), qc.blueF());

            Vector3 c;
            switch(shadeType) {
            case 1: { // Cosine like
                double s = lightMap.at(i, j);
                s *= s;
                c = 0.2 * Vector3(1.0, 1.0, 1.0) + 0.6 * s * cz + 0.2 * s * Vector3(1.0, 1.0, 1.0);
                break;
            }
            case 2: { // Normal
                Vector3 n = hf.Normal(i, j);
                double s = n * lightDirection; //Normalized(Vector3(-2.0, 1.0, 4.0));
                s = 0.5 * (1.0 + s);
                s *= s;
                s *= s;
                Vector3 cn = Math::Lerp(Vector3(52, 78, 118)/255.0, cz, s);
                c = 0.15 * Vector3(0.975, 0.975, 0.975) + 0.85 * cn;
                break;
            }
            case 3: { // Relief
                Vector3 n = hf.Normal(i, j);
                double s = n * lightDirection; //Normalized(Vector3(1.0, 0.5, 2.5));
                s = 0.5 * (1.0 + s);
                s *= s;
                s *= s;
                double t = Vector2(n) * Normalized(Vector2(1, 1));
                t = 0.5 * (1.0 + t);
                Vector3 cn = s * Math::Lerp(Vector3(0.65, 0.75, 0.85), Vector3(1.00, 0.95, 0.80), t);
                c = 0.25 * Vector3(0.975, 0.975, 0.975) + 0.25 * cz + 0.50 * cn;
                break;
            }
            default:
                c = cz;
                break;
            }

            if (shadows) {
                c = c * (0.7*shadowMap.at(i, j) + 0.3);
            }

            shading.setPixelColor(i, j, toQColor(c).rgb());
        }
    }

    return shading;
}


QImage singleLandformImage(const IntField2& lands, int v, const QColor& colorLand = Qt::red, const QColor& colorBg = Qt::white)
{
    QImage img(lands.getSizeX(), lands.getSizeY(), QImage::Format_RGBA8888);
    for (int i = 0; i < img.width(); i++) {
        for (int j = 0; j < img.height(); j++) {
            if (lands.at(i, j) == v) {
                img.setPixelColor(i, j, colorLand);
            }
            else {
                img.setPixelColor(i, j, colorBg);
            }
         }
    }
    return img;
}

QImage thresholdImage(const ScalarField2& field, int v, const QColor& colorPos = Qt::red, const QColor& colorNeg = Qt::white)
{
    QImage img(field.getSizeX(), field.getSizeY(), QImage::Format_RGBA8888);
    for (int i = 0; i < img.width(); i++) {
        for (int j = 0; j < img.height(); j++) {
            if (field.at(i, j) > v) {
                img.setPixelColor(i, j, colorPos);
            }
            else {
                img.setPixelColor(i, j, colorNeg);
            }
         }
    }
    return img;
}


ColorPalette MainWindow::getPalette(const ColorPalette& preferred) const
{
    if (!ui->fixed_palette->isChecked())
        return preferred;

    switch (ui->cb_palette->currentIndex()) {
        case 0: return ColorPalette::CoolWarm();
        case 1: return ColorPalette::Reds();
        case 2: return ColorPalette::Blues();
        default: return ColorPalette::CoolWarm();
    }
}

void getPaletteRange(const std::vector<double>& range, bool centered, double& rMin, double& rMax)
{
    if (centered) {
        double r = std::max(std::abs(range[0]), std::abs(range[1]));
        rMin = -r;
        rMax =  r;
    }
    else {
        rMin = range[0];
        rMax = range[1];
    }
}


void MainWindow::computeMetric()
{
    // LOCAL
    if (ui->dem_slopeGradient->isChecked()) {
        currentMetric = hf.GradientNorm();
        paletteDefault = ColorPalette::Reds();
        getPaletteRange(currentMetric.percentiles({0.001, 0.999}), false, paletteMin, paletteMax);
    }
    else if (ui->dem_slopeAvg->isChecked()) {
        currentMetric = hf.SlopeAverage();
        paletteDefault = ColorPalette::Reds();
        getPaletteRange(currentMetric.percentiles({0.001, 0.999}), false, paletteMin, paletteMax);
    }
    else if (ui->dem_aspect->isChecked()) {
        currentMetric = hf.Aspect();
        paletteDefault = ColorPalette::CoolWarm();
        paletteMin = 0;
        paletteMax = 2*Math::Pi;
    }
    else if (ui->dem_aspect_east->isChecked()) {
        currentMetric = hf.Aspect();
        for (int i = 0; i < currentMetric.getNumElements(); i++)
            currentMetric[i] = std::sin(currentMetric[i]);
        paletteDefault = ColorPalette::Reds();
        paletteMin = -1;
        paletteMax = 1;
    }
    else if (ui->dem_aspect_north->isChecked()) {
        currentMetric = hf.Aspect();
        for (int i = 0; i < currentMetric.getNumElements(); i++)
            currentMetric[i] = std::cos(currentMetric[i]);
        paletteDefault = ColorPalette::Reds();
        paletteMin = -1;
        paletteMax = 1;
    }

    // LAPLACIAN
    else if (ui->dem_laplacian->isChecked()) {
        currentMetric = hf.Laplacian();
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.01, 0.99}), true, paletteMin, paletteMax);
    }
    else if (ui->dem_fractLaplacian->isChecked()) {
        currentMetric = hf.FractLaplacian(ui->dem_fractLaplacian_s->value(), ui->dem_fractLaplacian_n->value());
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.01, 0.99}), true, paletteMin, paletteMax);
    }

    // ASYMMETRY
    else if (ui->dem_asymmetry->isChecked()) {
        currentMetric = hf.HillslopeAsymmetry(ui->dem_asymmetryW->value(), ui->dem_asymmetryDir->value(), ui->dem_asymmetryDirTol->value());
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.01, 0.99}), true, paletteMin, paletteMax);
    }

    // CURVATURES
    else if (ui->dem_curvMin->isChecked()) {
        currentMetric = hf.Curvature(HeightField::CurvatureType::MIN,
                                     ui->dem_curv_fromQuadrics->isChecked() ? ui->dem_curv_fromQuadrics_w->value() : 0);
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.1, 0.9}), true, paletteMin, paletteMax);
    }
    else if (ui->dem_curvMax->isChecked()) {
        currentMetric = hf.Curvature(HeightField::CurvatureType::MAX,
                                     ui->dem_curv_fromQuadrics->isChecked() ? ui->dem_curv_fromQuadrics_w->value() : 0);
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.1, 0.9}), true, paletteMin, paletteMax);
    }
    else if (ui->dem_curvMean->isChecked()) {
        currentMetric = hf.Curvature(HeightField::CurvatureType::MEAN,
                                     ui->dem_curv_fromQuadrics->isChecked() ? ui->dem_curv_fromQuadrics_w->value() : 0);
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.1, 0.9}), true, paletteMin, paletteMax);
    }
    else if (ui->dem_curvGaussian->isChecked()) {
        currentMetric = hf.Curvature(HeightField::CurvatureType::GAUSSIAN,
                                     ui->dem_curv_fromQuadrics->isChecked() ? ui->dem_curv_fromQuadrics_w->value() : 0);
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.01, 0.99}), true, paletteMin, paletteMax);
    }
    else if (ui->dem_curvProfile->isChecked()) {
        currentMetric = hf.Curvature(HeightField::CurvatureType::PROFILE,
                                     ui->dem_curv_fromQuadrics->isChecked() ? ui->dem_curv_fromQuadrics_w->value() : 0);
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.1, 0.9}), true, paletteMin, paletteMax);
    }
    else if (ui->dem_curvContour->isChecked()) {
        currentMetric = hf.Curvature(HeightField::CurvatureType::CONTOUR,
                                     ui->dem_curv_fromQuadrics->isChecked() ? ui->dem_curv_fromQuadrics_w->value() : 0);
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.1, 0.9}), true, paletteMin, paletteMax);
    }
    else if (ui->dem_curvTangent->isChecked()) {
        currentMetric = hf.Curvature(HeightField::CurvatureType::TANGENTIAL,
                                         ui->dem_curv_fromQuadrics->isChecked() ? ui->dem_curv_fromQuadrics_w->value() : 0);
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.1, 0.9}), true, paletteMin, paletteMax);
    }

    // RELIEF
    else if (ui->dem_localRelief->isChecked()) {
        currentMetric = hf.LocalRelief(ui->dem_localRelief_w->value());
        paletteDefault = ColorPalette::Reds();
        getPaletteRange(currentMetric.percentiles({0.001, 0.999}), false, paletteMin, paletteMax);
    }
    else if (ui->dem_localVariance->isChecked()) {
        currentMetric = hf.LocalVariance(ui->dem_localRelief_w->value());
        paletteDefault = ColorPalette::Reds();
        getPaletteRange(currentMetric.percentiles({0.001, 0.999}), false, paletteMin, paletteMax);
    }
    else if (ui->dem_areaRatio->isChecked()) {
        currentMetric = hf.AreaRatio(ui->dem_localRelief_w->value());
        paletteDefault = ColorPalette::Reds();
        paletteMin = 1;
        paletteMax = currentMetric.percentile(0.999);
    }
    else if (ui->dem_tpi->isChecked()) {
        currentMetric = hf.TopographicPositionIndexSAT(ui->dem_localRelief_w->value());
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.001, 0.999}), true, paletteMin, paletteMax);
    }
    else if (ui->dem_ruggedness->isChecked()) {
        currentMetric = hf.RuggednessIndex(ui->dem_localRelief_w->value());
        paletteDefault = ColorPalette::Reds();
        getPaletteRange(currentMetric.percentiles({0.001, 0.999}), false, paletteMin, paletteMax);
    }
    else if (ui->dem_surfaceRoughness->isChecked()) {
        currentMetric = hf.SurfaceRoughnessSAT(ui->dem_localRelief_w->value(), false);
        paletteDefault = ColorPalette::Reds();
        paletteMin = 0;
        paletteMax = currentMetric.percentile(0.999);
    }
    else if (ui->dem_surfaceRoughnessStd->isChecked()) {
        currentMetric = hf.SurfaceRoughnessSAT(ui->dem_localRelief_w->value(), true);
        paletteDefault = ColorPalette::Reds();
        paletteMin = 0;
        paletteMax = currentMetric.percentile(0.999);
    }

    // VISIBILITY
    else if (ui->dem_viewshedTotalIn->isChecked()) {
        double s = ui->viewshed_rescale->value();
        HeightField hfLow = hf.setResolution(int(s*hf.getSizeX()), int(s*hf.getSizeY()));
        ScalarField2 view;
        if (ui->dem_viewshedTotal_sample->isChecked() && ui->dem_viewshedTotal_numSamples->value() < hf.getNumElements())
            view = hfLow.ViewshedTotalSampled(ui->dem_viewshedTotal_numSamples->value(), false, ui->viewshed_offset->value());
        else
            view = hfLow.ViewshedTotal(false, ui->viewshed_offset->value());
        view *= 1.0/view.getNumElements();
        currentMetric = view.setResolution(hf.getSizeX(), hf.getSizeY());
        paletteDefault = ColorPalette::Reds();
        currentMetric.getRange(paletteMin, paletteMax);
        paletteMin = 0;
    }
    else if (ui->dem_viewshedTotalOut->isChecked()) {
        double s = ui->viewshed_rescale->value();
        HeightField hfLow = hf.setResolution(int(s*hf.getSizeX()), int(s*hf.getSizeY()));
        ScalarField2 view;
        if (ui->dem_viewshedTotal_sample->isChecked() && ui->dem_viewshedTotal_numSamples->value() < hf.getNumElements())
            view = hfLow.ViewshedTotalSampled(ui->dem_viewshedTotal_numSamples->value(), true, ui->viewshed_offset->value());
        else
            view = hfLow.ViewshedTotal(true, ui->viewshed_offset->value());
        view *= 1.0/view.getNumElements();
        currentMetric = view.setResolution(hf.getSizeX(), hf.getSizeY());
        paletteDefault = ColorPalette::Reds();
        currentMetric.getRange(paletteMin, paletteMax);
        paletteMin = 0;
    }
    else if (ui->dem_accessibility->isChecked()) {
        double r = ui->dem_accessibilityRadius->value();
        if (r <= 0) r = 2*hf.getDomain().radius();
        currentMetric = hf.Accessibility(r, ui->dem_accessibilitySamples->value());
        paletteDefault = ColorPalette::Reds();
        paletteMin = 0.5;
        paletteMax = 1.0;
    }
    else if (ui->dem_skyview->isChecked()) {
        double r = ui->dem_accessibilityRadius->value();
        if (r <= 0) r = 2*hf.getDomain().radius();
        currentMetric = hf.Accessibility(r, ui->dem_accessibilitySamples->value(), true);
        paletteDefault = ColorPalette::Reds();
        paletteMin = 0.5;
        paletteMax = 1.0;
    }
    else if (ui->dem_skyviewApprox->isChecked()) {
        currentMetric = hf.SkyViewFactorApproximation();
        paletteDefault = ColorPalette::Reds();
        paletteMin = 0.5;
        paletteMax = 1.0;
    }
    else if (ui->dem_posOpenness->isChecked()) {
        currentMetric = hf.Openness(true, ui->dem_opennessDist->value(), ui->dem_opennessDirs->value());
        paletteDefault = ColorPalette::CoolWarm();
        paletteMin = Math::Pi / 3;
        paletteMax = 2 * Math::Pi / 3;
    }
    else if (ui->dem_negOpenness->isChecked()) {
        currentMetric = hf.Openness(false, ui->dem_opennessDist->value(), ui->dem_opennessDirs->value());
        paletteDefault = ColorPalette::CoolWarm();
        paletteMin = Math::Pi / 3;
        paletteMax = 2 * Math::Pi / 3;
    }
    else if (ui->dem_sunlight->isChecked()) {
        currentMetric = hf.Sun(ui->dem_sunlightLatitude->value(), 14, 24, 7, 12);
        paletteDefault = ColorPalette::Reds();
        currentMetric.getRange(paletteMin, paletteMax);
        paletteMin = 0;
    }
    else if (ui->dem_dahi->isChecked()) {
        currentMetric = hf.DiurnalAnisotropicHeatIndex(Math::DegreeToRadian(ui->dem_dahiMax->value()));
        paletteDefault = ColorPalette::Reds();
        paletteMin = -0.65;
        paletteMax =  0.65;
    }

    // HYDROLOGY
    else if (ui->dem_streamArea->isChecked() || ui->dem_streamAreaLog->isChecked()
             || ui->dem_wetness->isChecked() || ui->dem_streamPower->isChecked())
    {
        std::vector<Index2> sources;
        if (ui->dem_streamArea_uniqueSource->isChecked()) {
            sources.push_back(Index2(ui->dem_streamArea_source_x->value(), ui->dem_streamArea_source_y->value()));
        }

        ScalarField2 sa;
        switch (ui->dem_streamArea_flowAlgorithm->currentIndex()) {
            case 0: sa = hf.StreamAreaD8(sources); break;
            case 1: sa = hf.StreamAreaMFD(sources, ui->dem_streamArea_mfdExp->value()); break;
            case 2: sa = hf.StreamAreaDinf(sources); break;
            case 3: sa = hf.StreamAreaRho8(sources, ui->dem_streamArea_rho8iters->value()); break;
            case 4: sa = hf.StreamAreaKinematic(sources); break;
        }

        if (ui->dem_wetness->isChecked()) {
            currentMetric = hf.WetnessIndex(sa);
            getPaletteRange(currentMetric.percentiles({0.01, 0.99}), false, paletteMin, paletteMax);
        }
        else if (ui->dem_streamPower->isChecked()) {
            currentMetric = hf.StreamPower(sa, ui->dem_streamPower_m->value(), ui->dem_streamPower_n->value());
            getPaletteRange(currentMetric.percentiles({0.01, 0.99}), false, paletteMin, paletteMax);
        }
        else {
            for (int i = 0; i < sa.getNumElements(); i++) {
                if (ui->dem_streamArea->isChecked())
                    sa[i] = std::pow(sa.at(i), ui->dem_streamArea_exp->value());
                else if (ui->dem_streamAreaLog->isChecked())
                    sa[i] = std::log(sa.at(i));
            }
            currentMetric = sa;
            paletteMin = 1;
            paletteMax = currentMetric.percentile(0.999);
        }
        paletteDefault = ColorPalette::Blues();
    }
    else if (ui->dem_depthToWater->isChecked()) {
        terrainAnalysis->computeRivers(ui->rivers_cit_s->value(), ui->rivers_cit_t->value());
        currentMetric = terrainAnalysis->DepthToWater();
        paletteDefault = ColorPalette::Blues();
        paletteMin = 0;
        paletteMax = currentMetric.percentile(0.999);
    }
    else if (ui->dem_hand->isChecked()) {
        terrainAnalysis->computeRivers(ui->rivers_cit_s->value(), ui->rivers_cit_t->value());
        currentMetric = terrainAnalysis->HeightAboveNearestDrainage();
        paletteDefault = ColorPalette::Blues();
        paletteMin = 0;
        paletteMax = currentMetric.percentile(0.999);
    }
    else if (ui->dem_riverDistEuclidean->isChecked()) {
        terrainAnalysis->computeRivers(ui->rivers_cit_s->value(), ui->rivers_cit_t->value());
        currentMetric = terrainAnalysis->DistanceToNearestDrainageEuclidean();
        paletteDefault = ColorPalette::Blues();
        paletteMin = 0;
        paletteMax = currentMetric.percentile(0.999);
    }
    else if (ui->dem_riverDistFlow->isChecked()) {
        terrainAnalysis->computeRivers(ui->rivers_cit_s->value(), ui->rivers_cit_t->value());
        currentMetric = terrainAnalysis->DistanceToNearestDrainageFlow();
        paletteDefault = ColorPalette::Blues();
        paletteMin = 0;
        paletteMax = currentMetric.percentile(0.999);
    }
    else if (ui->dem_branchLength->isChecked()) {
        ScalarField2 bmax = terrainAnalysis->MaximumBranchLength();
        double e = ui->dem_branchLength_exp->value();
        for (int i = 0; i < currentMetric.getNumElements(); i++) {
            bmax[i] = std::pow(bmax.at(i), e);
        }
        currentMetric = bmax;
        paletteDefault = ColorPalette::Blues();
        paletteMin = 0;
        paletteMax = currentMetric.percentile(0.999);
    }

    else if (ui->dem_peakedness->isChecked()) {
        currentMetric = hf.PeakPercentage(ui->dem_orometryRadius->value()*1000);
        paletteDefault = ColorPalette::Reds();
        paletteMin = 0;
        paletteMax = 1;
    }
    else if (ui->dem_ors->isChecked()) {
        currentMetric = hf.ORS(ui->dem_orometryRadius->value()*1000);
        paletteDefault = ColorPalette::Reds();
        getPaletteRange(currentMetric.percentiles({0.001, 0.999}), false, paletteMin, paletteMax);
    }
    else if (ui->dem_jut->isChecked()) {
        if (ui->dem_jut_curvature->isChecked())
            currentMetric = hf.JutCurved(ui->dem_orometryRadius->value()*1000, ui->dem_jut_planetR->value()*1000);
        else
            currentMetric = hf.JutPlanar(ui->dem_orometryRadius->value()*1000);
        paletteDefault = ColorPalette::Reds();
        getPaletteRange(currentMetric.percentiles({0.001, 0.999}), false, paletteMin, paletteMax);
    }
    else if (ui->dem_rut->isChecked()) {
        if (ui->dem_jut_curvature->isChecked())
            currentMetric = hf.RutCurved(ui->dem_orometryRadius->value()*1000, ui->dem_jut_planetR->value()*1000);
        else
            currentMetric = hf.RutPlanar(ui->dem_orometryRadius->value()*1000);
        paletteDefault = ColorPalette::Reds();
        getPaletteRange(currentMetric.percentiles({0.001, 0.999}), false, paletteMin, paletteMax);
    }
    else if (ui->dem_zreduced_converging->isChecked()) {
        currentMetric = hf.AngleReducedHeight(ui->dem_jut_planetR->value()*1000,
                                              ui->dem_zreduced_x->value(), ui->dem_zreduced_y->value(), true);
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.001, 0.999}), true, paletteMin, paletteMax);
    }
    else if (ui->dem_zreduced_diverging->isChecked()) {
        currentMetric = hf.AngleReducedHeight(ui->dem_jut_planetR->value()*1000,
                                              ui->dem_zreduced_x->value(), ui->dem_zreduced_y->value(), false);
        paletteDefault = ColorPalette::CoolWarm();
        getPaletteRange(currentMetric.percentiles({0.001, 0.999}), true, paletteMin, paletteMax);
    }
    else if (ui->dem_deng2008_peakProto->isChecked() || ui->dem_deng2008_peakSimil->isChecked()) {
        ScalarField2 peakProto;
        ScalarField2 peakSimil = hf.Peakness(peakProto, ui->dem_deng2008_peakPercent->value()/100.0);
        if (ui->dem_deng2008_peakProto->isChecked())
            currentMetric = peakProto;
        else
            currentMetric = peakSimil;
        paletteDefault = ColorPalette::Reds();
        currentMetric.getRange(paletteMin, paletteMax);
    }

    else if (ui->dem_nni->isChecked()) {
        double terrainNNI;
        terrainAnalysis->computeDivideTree(ui->ridges_divtree_promCutoff->value());
        currentMetric = terrainAnalysis->NearestNeighborIndex(ui->dem_nni_radius->value()*1000, terrainNNI);

        if (ui->dem_nni_relative->isChecked()) {
            for (int i = 0; i < currentMetric.getNumElements(); i++) {
                currentMetric[i] -= terrainNNI;
            }
            paletteDefault = ColorPalette::CoolWarm();
            getPaletteRange(currentMetric.percentiles({0.001, 0.999}), true, paletteMin, paletteMax);
        }
        else {
            paletteDefault = ColorPalette::Reds();
            getPaletteRange(currentMetric.percentiles({0.001, 0.999}), false, paletteMin, paletteMax);
        }
        ui->labelTerrainNNI->setText(QString::number(terrainNNI));
    }
}

void MainWindow::computeBaseTexture(bool updateMetric)
{
    shadeType = 1;
    const QColor whiteFill = QColor(232, 232, 232);

    // first the descriptors with a fixed palette, then the general case with UI-configurable palette

    // plain model
    if (ui->dem_uniform->isChecked()) {
        baseTexture = QImage(hf.getSizeX(), hf.getSizeY(), QImage::Format_RGBA8888);
        baseTexture.fill(whiteFill);
        shadeType = 1;
        currentMetric = hf;
    }

    // hypsometric gradient
    else if (ui->dem_elevation->isChecked()) {
        double zmin = hf.getBox().getMin()[2];
        double zmax = hf.getBox().getMax()[2];
        baseTexture = hf.createImage(zmin, zmin + (zmax - zmin)/0.75, ColorPalette::Relief());
        shadeType = 3;
        currentMetric = hf;
    }

    // water level
    else if (ui->dem_waterlevel->isChecked()) {
        double waterLevel = ui->dem_waterlevel_m->value();
        bool relief = ui->dem_waterlevel_relief->isChecked();
        baseTexture = QImage(hf.getSizeX(), hf.getSizeY(), QImage::Format_RGBA8888);
        baseTexture.fill(whiteFill);
        for (int i = 0; i < baseTexture.width(); i++) {
            for (int j = 0; j < baseTexture.height(); j++) {
                if (hf.at(i, j) < waterLevel) {
                    if (relief) {
                        double u = std::max(0.0, 1 - (waterLevel - hf.at(i, j)) / 1000.0);
                        Vector3 col = Math::Lerp(Vector3(18, 53, 102), Vector3(52, 135, 255), u)/255.0;
                        baseTexture.setPixelColor(i, j, toQColor(col));
                    }
                    else {
                        baseTexture.setPixelColor(i, j, QColor(48, 128, 255, 255));
                    }
                }
            }
        }
        shadeType = 1;
        currentMetric = hf;
    }

    // normal map
    else if (ui->dem_normals->isChecked()) {
        baseTexture = QImage(hf.getSizeX(), hf.getSizeY(), QImage::Format_RGBA8888);
        for (int i = 0; i < hf.getSizeX(); i++) {
            for (int j = 0; j < hf.getSizeY(); j++) {
                Vector3 n = 0.5 * (hf.Normal(i, j) + Vector3(1));
                baseTexture.setPixelColor(i, j, toQColor(Vector3(n[0], n[1], n[2])));
            }
        }
        shadeType = 0;
        currentMetric = hf;
    }

    // binary viewshed
    else if (ui->dem_viewshedLoc->isChecked()) {
        double s = ui->viewshed_rescale->value();
        HeightField hfLow = hf.setResolution(int(s*hf.getSizeX()), int(s*hf.getSizeY()));

        int vpx = ui->viewshed_location_x->value();
        int vpy = ui->viewshed_location_y->value();
        double offset = ui->viewshed_offset->value();

        int visCells;
        ScalarField2 view = hfLow.Viewshed(hf.Vertex(vpx, vpy) + Vector3(0, 0, offset), visCells);

        baseTexture = QImage(hfLow.getSizeX(), hfLow.getSizeY(), QImage::Format_RGBA8888);
        baseTexture.fill(whiteFill);
        for (int i = 0; i < baseTexture.width(); i++) {
            for (int j = 0; j < baseTexture.height(); j++) {
                if (view.at(i, j) > 0) baseTexture.setPixelColor(i, j, Qt::red);
            }
        }
        baseTexture = baseTexture.scaled(QSize(hf.getSizeX(), hf.getSizeY()), Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
        baseTexture.setPixelColor(vpx, vpy, Qt::green);
        currentMetric = view.setResolution(hf.getSizeX(), hf.getSizeY());
    }

    // landforms
    else if (ui->dem_dikauWoodClass->isChecked()) {
        ui->dem_showSingleLandform_id->setMaximum(namesLandDW.size() - 1);
        ui->dem_landformName->setText(namesLandDW[ui->dem_showSingleLandform_id->value()]);

        IntField2 lands = hf.LandformsDikauWood(ui->dem_dikauWoodSolver_w->value(), 0.1, 0.0001);
        if (ui->dem_showSingleLandform->isChecked())
            baseTexture = singleLandformImage(lands, ui->dem_showSingleLandform_id->value());
        else
            baseTexture = lands.createImage(paletteLandDW);
        currentMetric = lands.toScalarField();
    }
    else if (ui->dem_fisher2004->isChecked() || ui->dem_fisher2004entropy->isChecked()) {
        ScalarField2 entropy(hf.getDomain(), hf.getSizeX(), hf.getSizeY());
        double tCurv = 0.0001;
        double tSlope = 0.1;
        int minScale = ui->dem_fisher2004_minScale->value();
        int maxScale = ui->dem_fisher2004_maxScale->value();
        IntField2 lands = hf.LandformsFuzzyDW(entropy, minScale, maxScale, tSlope, tCurv);

        if (ui->dem_fisher2004entropy->isChecked()) {
            ui->dem_showSingleLandform_id->setMaximum(0);
            ui->dem_landformName->setText("-");

            baseTexture = entropy.createImage(ColorPalette::Reds());
            currentMetric = entropy;
        }
        else {
            ui->dem_showSingleLandform_id->setMaximum(namesLandDW.size() - 1);
            ui->dem_landformName->setText(namesLandDW[ui->dem_showSingleLandform_id->value()]);

            if (ui->dem_showSingleLandform->isChecked())
                baseTexture = singleLandformImage(lands, ui->dem_showSingleLandform_id->value());
            else
                baseTexture = lands.createImage(paletteLandDW);
            currentMetric = lands.toScalarField();
        }
    }
    else if (ui->dem_tpiClass->isChecked()) {
        ui->dem_showSingleLandform_id->setMaximum(namesLandTPI.size() - 1);
        ui->dem_landformName->setText(namesLandTPI[ui->dem_showSingleLandform_id->value()]);

        IntField2 lands = hf.LandformsTPI(ui->dem_tpiClassSmall->value(), ui->dem_tpiClassLarge->value());
        if (ui->dem_showSingleLandform->isChecked())
            baseTexture = singleLandformImage(lands, ui->dem_showSingleLandform_id->value());
        else
            baseTexture = lands.createImage(paletteLandTPI);
        currentMetric = lands.toScalarField();
    }
    else if (ui->dem_geomorphons->isChecked()) {
        ui->dem_showSingleLandform_id->setMaximum(namesLandGeomorphons.size() - 1);
        ui->dem_landformName->setText(namesLandGeomorphons[ui->dem_showSingleLandform_id->value()]);

        IntField2 lands = hf.Geomorphons(ui->dem_geomorphons_km->value() * 1000.0);
        if (ui->dem_showSingleLandform->isChecked())
            baseTexture = singleLandformImage(lands, ui->dem_showSingleLandform_id->value());
        else
            baseTexture = lands.createImage(paletteLandGeomorphons);
        currentMetric = lands.toScalarField();
    }
    else if (ui->dem_tophat->isChecked()) {
        ui->dem_showSingleLandform_id->setMaximum(0);
        ui->dem_landformName->setText("-");

        int w = ui->dem_tophat_w->value();
        double tPeak = ui->dem_tophat_peak->value();
        double tValley = ui->dem_tophat_valley->value();
        currentMetric = hf.BlackWhiteTopHatTransform(w, tPeak, tValley);
        baseTexture = currentMetric.createImage(ColorPalette::CoolWarm());
    }

    // all the rest, the palette is controlled via UI so we do not want
    // to recompute the metric every time the range is modified
    else {
        if (updateMetric) {
            computeMetric();
            updatePaletteRangeSelectors();
            updateCursorInfoTexts();
        }
        baseTexture = currentMetric.createImage(ui->sb_paletteMin->value() * paletteScale,
                                                ui->sb_paletteMax->value() * paletteScale,
                                                getPalette(paletteDefault));

        // some metrics still require a special treatment
        if (ui->dem_streamPower->isChecked() && ui->dem_streamPower_threshold->isChecked()) {
            baseTexture = thresholdImage(currentMetric, ui->dem_streamPower_t->value());
        }
        else if (ui->dem_zreduced_converging->isChecked() || ui->dem_zreduced_diverging->isChecked()) {
            baseTexture.setPixelColor(ui->dem_zreduced_x->value(), ui->dem_zreduced_y->value(), Qt::green);
        }
    }

}



inline QPointF getScaledPoint(const Offsets& off, int scale = 1)
{
  return QPointF(scale * (off.x() + 0.5), scale * (off.y() + 0.5));
}

inline QPointF getScaledPoint(const Index2& off, int scale = 1)
{
  return QPointF(scale * (off.x() + 0.5), scale * (off.y() + 0.5));
}

QColor interpolateColor(const QColor &color1, const QColor &color2, double t) {
    t = Math::Clamp(t, 0.0, 1.0);
    int r = static_cast<int>(color1.red() + (color2.red() - color1.red()) * t);
    int g = static_cast<int>(color1.green() + (color2.green() - color1.green()) * t);
    int b = static_cast<int>(color1.blue() + (color2.blue() - color1.blue()) * t);
    return QColor(r, g, b);
}


void MainWindow::computeRiversLayer()
{
    riversLayer = QImage(hf.getSizeX(), hf.getSizeY(), QImage::Format_RGBA8888);
    riversLayer.fill(QColor(0, 0, 0, 0));

    double cit_s = ui->rivers_cit_s->value();
    double cit_t = ui->rivers_cit_t->value();
    std::vector<RiverTree*> rivers = terrainAnalysis->computeRivers(cit_s, cit_t);
    TerrainFlowD8* flowD8 = terrainAnalysis->getFlowDirs();

    ScalarField2 riverStreamLength(hf.getDomain(), hf.getSizeX(), hf.getSizeY(), 0);
    ScalarField2 riverStreamArea(hf.getDomain(), hf.getSizeX(), hf.getSizeY(), 0);

    // draw basins
    if (ui->showBasins->isChecked()) {

        int basinsAlpha = 192;
        std::vector<QColor> basinsPalette = {
            QColor(198,219,239,basinsAlpha),
            QColor(158,202,225,basinsAlpha),
            QColor(107,174,214,basinsAlpha),
            QColor(49,130,189,basinsAlpha)
        };
        QColor basinBlue(158,202,225,basinsAlpha);

        // single basin from coords
        if (ui->rivers_drainageBasinCoords->isChecked()) {
            Index2 pSelect(ui->rivers_basin_x->value(), ui->rivers_basin_y->value());
            Index2 p = pSelect;
            if (ui->rivers_basin_useRiver->isChecked()) {
                const IntField2& mask = terrainAnalysis->RiversMask();
                while (mask.isValidCell(p) && mask.at(p) == 0) {
                    if (!flowD8->hasDownstream(p)) break;
                    p = flowD8->downstream(p);
                }
            }

            DrainageBasin basin = flowD8->getDrainageBasin(p);
            for (int i = 0; i < hf.getSizeX(); i++) {
                for (int j = 0; j < hf.getSizeY(); j++) {
                    if (basin.inBasin(i, j)) riversLayer.setPixelColor(i, j, basinBlue);
                }
            }
            riversLayer.setPixelColor(pSelect.x(), pSelect.y(), Qt::green);
        }
        // multiple basins
        else {

            // get basin endpoints
            std::vector<Index2> basinPoints;
            if (ui->rivers_drainageBasinAll->isChecked()) {
                for (const RiverTree* r : rivers) {
                    basinPoints.push_back(r->getCell(0));
                    r->markRiverLength(riverStreamLength);
                    r->markRiverStreamArea(riverStreamArea);
                }
            }
            else if (ui->rivers_drainageBasinStrahler->isChecked()) {
                int strahler = ui->rivers_drainageBasinStrahler_order->value();
                std::queue<const RiverTree*> Q;
                for (const RiverTree* r : rivers) {
                    if (r->getStrahlerIndex() >= strahler)
                        Q.push(r);
                }
                while (!Q.empty()) {
                    const RiverTree* r = Q.front();
                    Q.pop();
                    if (r->getStrahlerIndex() > strahler) {
                        for (const RiverTree* p : r->getParentRivers())
                            Q.push(p);
                    }
                    else if (r->getStrahlerIndex() == strahler) {
                        if (r->getNumCells() > 1) {
                            // if we have at least 2 cells, we might come from a junction (stored at 0)
                            // and we want to omit this cell, since it represents already strahler+1
                            basinPoints.push_back(r->getCell(1));
                        }
                        else {
                            basinPoints.push_back(r->getCell(0));
                        }
                        r->markRiverLength(riverStreamLength);
                        r->markRiverStreamArea(riverStreamArea);
                    }
                }
            }

            // get basins
            std::vector<DrainageBasin> basins;
            std::vector<double> basinMetric;

            for (unsigned int i = 0; i < basinPoints.size(); i++) {
                Index2 p = basinPoints[i];
                DrainageBasin basin = flowD8->getDrainageBasin(p);
                basin.computeBasinMetrics(hf, riverStreamLength, riverStreamArea);
                basins.push_back(basin);

                double val;
                switch (ui->rivers_drainageBasin_color->currentIndex()) {
                case 0: // Id
                    val = i;
                    break;
                case 1: // Area
                    val = basin.getArea();
                    break;
                case 2: // Perimeter
                    val = basin.getPerimeter();
                    break;
                case 3: // Circularity
                    val = basin.getCircularity();
                    break;
                case 4: // Compactness
                    val = basin.getCompactness();
                    break;
                case 5: // MRN
                    val = basin.getMeltonRuggednessIndex();
                    break;
                case 6: // Hypsommetric integral
                    val = basin.getHypsommetricIntegral();
                    break;
                case 7: // Drainage area
                    val = basin.getDrainageDensity();
                    break;
                case 8: // Relief ratio
                    val = basin.getReliefRatio();
                    break;
                case 9: // Form factor
                    val = basin.getFormFactor();
                    break;
                default: // unknown
                    break;
                }

                basinMetric.push_back(val);
            }

            bool perim = ui->rivers_drainageBasin_perimeter->isChecked();
            for (unsigned int i = 0; i < basins.size(); i++) {
                const DrainageBasin& basin = basins[i];
                QColor color;
                if (ui->rivers_drainageBasin_color->currentIndex() == 0) {
                    color = basinsPalette[i%basinsPalette.size()];
                } else {                    
                    std::vector<double> sortedVals(basinMetric);
                    std::sort(sortedVals.begin(), sortedVals.end());
                    double vmin = sortedVals[int(0.1*sortedVals.size())];
                    double vmax = sortedVals[int(0.9*sortedVals.size())];
                    double t = (vmax - vmin) > 0 ? Math::Clamp((basinMetric[i] - vmin)/(vmax - vmin), 0.0, 1.0) : 0;
                    color = toQColor(ColorPalette::CoolWarm().getColor(t));
                    color.setAlpha(basinsAlpha);
                }
                for (int i = 0; i < hf.getSizeX(); i++) {
                    for (int j = 0; j < hf.getSizeY(); j++) {
                        if (perim && basin.inPerimeter(i, j))
                            riversLayer.setPixelColor(i, j, QColor(64,64,64,basinsAlpha));
                        else if (basin.inBasin(i, j))
                            riversLayer.setPixelColor(i, j, color);
                    }
                }
            }
        }
    }

    riversLayer = riversLayer.scaled(hf.getSizeX()*scaleLayerRivers, hf.getSizeY()*scaleLayerRivers,
                                     Qt::KeepAspectRatio, Qt::FastTransformation);

    if (!ui->showRivers->isChecked()) return;


    QPainter painter(&riversLayer);
    painter.setPen(QPen(Qt::blue, 1.0));
    double riverNodeSize = 0.75*scaleLayerRivers;
    QColor riverColor1(49, 201, 227);
    QColor riverColor2(78, 78, 192);

    bool simplePath = ui->showRiversSimple->isChecked();
    bool drawSources = ui->showRiversSources->isChecked();
    bool drawJunctions = ui->showRiversJunctions->isChecked();
    bool stopAtSea = ui->rivers_stopSea->isChecked();
    double seaLevel = ui->rivers_stopSea_elev->value();

    std::vector<QPointF> sources;
    std::vector<QPointF> junctions;
    std::vector<std::vector<QPointF> > segmentPaths;
    std::vector<RiverTree*> segmentRiver;
    double maxStreamArea = 0;
    double maxSinuosity = 0;
    double maxMeanSlope = 0;
    int maxStrahler = 1;

    for (RiverTree* r : rivers) {
        std::queue<RiverTree*> Q;
        Q.push(r);
        while (!Q.empty()) {
            RiverTree* river = Q.front();
            Q.pop();

            std::vector<QPointF> path;
            if (simplePath) {
                path.push_back(getScaledPoint(river->getStreamCells().front(), scaleLayerRivers));
                path.push_back(getScaledPoint(river->getStreamCells().back(), scaleLayerRivers));
            }
            else {
                path.resize(river->getNumCells());
                for (int i = 0; i < river->getNumCells(); i++) {
                    path[i] = getScaledPoint(river->getCell(i), scaleLayerRivers);
                }
            }
            segmentPaths.push_back(path);
            segmentRiver.push_back(river);
            maxStreamArea = std::max(maxStreamArea, river->getStreamArea().back());
            maxSinuosity  = std::max(maxSinuosity,  river->getSinuosity());
            maxMeanSlope  = std::max(maxMeanSlope,  Math::Mean(river->getStreamSlope()));
            maxStrahler   = std::max(maxStrahler,   river->getStrahlerIndex());

            // record source/junction circles
            if (river->getNumParentRivers() > 0) {
                junctions.push_back(path.back());
            }
            if (river->getNumParentRivers() == 0){
                sources.push_back(path.back());
            }

            // process parent river segments
            for (RiverTree* parent : river->getParentRivers()) {
                Q.push(parent);
            }
        }
    }


    // draw paths
    for (unsigned int i = 0; i < segmentPaths.size(); i++) {

        const std::vector<QPointF>& path = segmentPaths[i];
        const RiverTree* river = segmentRiver[i];
        QColor pathColor = riverColor2;
        double pathWidthDefault = 0.5 * scaleLayerRivers * ui->rivers_scale->value();
        double pathWidth = pathWidthDefault;

        int renderWidthType = ui->rivers_width->currentIndex();
        int renderColorType = ui->rivers_color->currentIndex();

        switch (renderWidthType) {
            case 0: // uniform
                pathWidth *= 1.5;
                break;
            case 1: // stream area
                pathWidth *= std::max(0.5, 6.0 * sqrt(river->getStreamArea().back()) / sqrt(maxStreamArea));
                break;
            case 2: // strahler
                pathWidth *= double(1 << (river->getStrahlerIndex() - 1));
                break;
            default: break;
        }
        double ct;
        switch (renderColorType) {
            case 0: break;
            case 1: // stream area
                ct = sqrt(river->getStreamArea().back()) / sqrt(maxStreamArea);
                pathColor = interpolateColor(riverColor1, riverColor2, ct);
                break;
            case 2: // strahler
                ct = double(river->getStrahlerIndex() - 1) / double(maxStrahler - 1);
                pathColor = interpolateColor(riverColor1, riverColor2, ct);
                break;
            case 3: // avg slope
                ct = Math::Mean(river->getStreamSlope()) / maxMeanSlope;
                pathColor = interpolateColor(riverColor1, riverColor2, ct);
                break;
            case 4: // sinuosity
                ct = (river->getSinuosity() - 1) / (maxSinuosity - 1);
                pathColor = interpolateColor(riverColor1, riverColor2, ct);
                break;
            default: break;
        }

        if (ui->showSegmentSmooth->isChecked()) {
            for (int i = 0; i < int(path.size()) - 1; ++i) {

                if (stopAtSea && hf.at(river->getCell(i)) < seaLevel)
                    continue;

                if (renderWidthType == 1) {
                    pathWidth = pathWidthDefault * std::max(0.5, 6.0 * sqrt(river->getStreamArea()[i]) / sqrt(maxStreamArea));
                }

                if (renderColorType == 1) {
                    ct = sqrt(river->getStreamArea()[i]) / sqrt(maxStreamArea);
                    pathColor = interpolateColor(riverColor1, riverColor2, ct);
                }
                else if (renderColorType == 3) {
                    ct = river->getStreamSlope()[i] / maxMeanSlope;
                    pathColor = interpolateColor(riverColor1, riverColor2, ct);
                }

                painter.setPen(QPen(QBrush(pathColor), pathWidth, Qt::SolidLine, Qt::RoundCap));
                painter.drawLine(path[i], path[i + 1]);
            }
        }
        else {
            unsigned int i0 = 0;
            if (stopAtSea) {
                while (i0 < path.size() && hf.at(river->getCell(i0)) < seaLevel) i0++;
            }
            if (i0 < path.size() - 1) {
                painter.setPen(QPen(QBrush(pathColor), pathWidth, Qt::SolidLine, Qt::RoundCap));
                painter.drawPolyline(&path[i0], int(path.size() - i0));
            }
        }
    }

    // draw source/junction circles
    if (drawSources) {
        for (const QPointF& p : sources) {
            painter.setPen(QPen(QBrush(QColor(0, 0, 0)), 2));
            painter.setBrush(QBrush(QColor(160, 228, 255)));
            painter.drawEllipse(p, riverNodeSize, riverNodeSize);
        }
    }
    if (drawJunctions) {
        for (const QPointF& p : junctions) {
            painter.setPen(QPen(QBrush(QColor(0, 0, 0)), 2));
            painter.setBrush(QBrush(QColor(64, 78, 204)));
            painter.drawEllipse(p, riverNodeSize, riverNodeSize);
        }
    }
}


void MainWindow::computeRidgesLayer()
{
    ridgesLayer = QImage(hf.getSizeX() * scaleLayerRidges, hf.getSizeY() * scaleLayerRidges, QImage::Format_RGBA8888);
    ridgesLayer.fill(QColor(0, 0, 0, 0));

    if (!ui->showRidges->isChecked()) return;

    int ridgeWidth = scaleLayerRidges;
    int peakSize = 8;
    int colSize = 4;

    QPainter painter(&ridgesLayer);
    painter.setPen(QPen(Qt::red, ridgeWidth));

    if (ui->ridges_drainage->isChecked()) {
        const double accThres = ui->ridges_drainageThreshold->value();
        ScalarField2 acc = hf.StreamAreaD8();
        QImage ridgeImg(acc.getSizeX(), acc.getSizeY(), QImage::Format_ARGB32);
        for (int i = 0; i < acc.getSizeX(); i++) {
            for (int j = 0; j < acc.getSizeY(); j++) {
                if (acc.at(i, j) <= accThres)
                    ridgeImg.setPixelColor(i, j, QColor(255, 0, 0, 255));
                else
                    ridgeImg.setPixelColor(i, j, QColor(0, 0, 0, 0));
            }
        }
        painter.drawImage(QPoint(0, 0), ridgeImg.scaled(ridgesLayer.size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));
    }
    else if (ui->ridges_laplacian->isChecked()) {
        const double lapThres = ui->ridges_laplacianThreshold->value();
        ScalarField2 lap = hf.Laplacian();
        QImage ridgeImg(lap.getSizeX(), lap.getSizeY(), QImage::Format_ARGB32);
        for (int i = 0; i < lap.getSizeX(); i++) {
            for (int j = 0; j < lap.getSizeY(); j++) {
                if (lap.at(i, j) <= lapThres)
                    ridgeImg.setPixelColor(i, j, QColor(255, 0, 0, 255));
                else
                    ridgeImg.setPixelColor(i, j, QColor(0, 0, 0, 0));
          }
        }
        painter.drawImage(QPoint(0, 0), ridgeImg.scaled(ridgesLayer.size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));
    }
    else if (ui->ridges_curvature->isChecked()) {
        const double curvThres = ui->ridges_curvatureThreshold->value();
        ScalarField2 curv = hf.Curvature(HeightField::CurvatureType::MIN, 0);
        QImage ridgeImg(curv.getSizeX(), curv.getSizeY(),  QImage::Format_ARGB32);
        for (int i = 0; i < curv.getSizeX(); i++) {
            for (int j = 0; j < curv.getSizeY(); j++) {
                if (curv.at(i, j) <= curvThres)
                    ridgeImg.setPixelColor(i, j, QColor(255, 0, 0, 255));
                else
                    ridgeImg.setPixelColor(i, j, QColor(0, 0, 0, 0));
            }
        }
        painter.drawImage(QPoint(0, 0), ridgeImg.scaled(ridgesLayer.size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));
    }

    else if (ui->ridges_PPA->isChecked()) {
        const PPA* ridgesPPA = terrainAnalysis->computePPA(ui->ridges_PPA_profileLength->value());
        const std::vector<PPA::RidgeSegment>& ridges = ridgesPPA->getSegmentsInPrunedMST();
        for (const PPA::RidgeSegment& s : ridges) {
            QPointF p1 = getScaledPoint(s.first, scaleLayerRidges);
            QPointF p2 = getScaledPoint(s.second, scaleLayerRidges);
            painter.drawLine(p1, p2);
        }
        if (ui->ridges_showPPAdebug->isChecked()) {
            painter.setPen(Qt::NoPen);
            // ridge candidate points
            painter.setBrush(QBrush(QColor(128, 0, 0)));
            for (const Index2& v : ridgesPPA->getRidgeCandidates()) {
                painter.drawEllipse(getScaledPoint(v, scaleLayerRidges), 3, 3);
            }
            // ridge connections
            painter.setPen(QPen(QColor(208, 128, 0), 1));
            for (const PPA::RidgeSegment& s : ridgesPPA->getSegments()) {
                QPointF p1 = getScaledPoint(s.first, scaleLayerRidges);
                QPointF p2 = getScaledPoint(s.second, scaleLayerRidges);
                painter.drawLine(p1, p2);
            }
            // reliable ridges
            painter.setPen(QPen(QColor(232, 208, 0), 1));
            for (const PPA::RidgeSegment& s : ridgesPPA->getReliableSegments()) {
                QPointF p1 = getScaledPoint(s.first, scaleLayerRidges);
                QPointF p2 = getScaledPoint(s.second, scaleLayerRidges);
                painter.drawLine(p1, p2);
            }
        }
    }

    else if (ui->ridges_divtree->isChecked()) {

        DivideTree* prunedDivTree = terrainAnalysis->computeDivideTree(ui->ridges_divtree_promCutoff->value());
        ui->ridges_showIsol_peak->setMaximum(int(prunedDivTree->peaks().size())-1);
        ui->labelNumPeaks->setText(QString::number(prunedDivTree->peaks().size()) + " peaks");
        bool ignoreSea = ui->ridges_ignoreSea->isChecked();
        double seaLevel = ui->ridges_ignoreSea_elev->value();

        std::vector<int> ridgeLevel;
        int maxLevel = 0;
        if (ui->ridges_level->isChecked()) {
            if (ui->ridges_level_type->currentIndex() == 0) {       // path length order
                terrainAnalysis->computeRidgesLevelRak2019(ridgeLevel, maxLevel);
            }
            else if (ui->ridges_level_type->currentIndex() == 1) {  // topological order
                terrainAnalysis->computeRidgesLevelScherler2020(ridgeLevel, maxLevel);
            }
        }

        std::vector<std::vector<Index2>> ridgeLines;
        std::vector<std::pair<int, int>> ridgePeaks;
        std::vector<int> ridgeSaddles;
        terrainAnalysis->getDivideTreeRidges(ui->ridges_divtree_detailed->isChecked(), ridgeLines, ridgePeaks, ridgeSaddles);

        double linew = ui->ridges_level_scale->value();

        for (int i = 0; linew > 0 && i < int(ridgeLines.size()); i++) {

            double sElev = prunedDivTree->getSaddle(ridgeSaddles[i]).elevation;
            if (ignoreSea && sElev <= seaLevel) continue;

            if (ui->ridges_level->isChecked()) {
                // path length order
                if (ui->ridges_level_type->currentIndex() == 0) {
                    int ridgeScale = maxLevel + 1 - ridgeLevel[i];
                    painter.setPen(QPen(QColor(int(255 * ((ridgeScale - 1) / float(maxLevel - 1))), 0, 0),
                                        linew * ridgeScale));
                }
                // topological order
                else if (ui->ridges_level_type->currentIndex() == 1) {
                    double ridgeScale = std::log(double(ridgeLevel[i])) + 1;
                    double ridgeLogMax = std::log(double(maxLevel)) + 1;
                    painter.setPen(QPen(QColor(int(255 * ridgeScale / ridgeLogMax), 0, 0),
                                        linew * ridgeScale));
                }
            }
            else {
                painter.setPen(QPen(Qt::red, 0.5 * linew * ridgeWidth));
            }

            const std::vector<Index2>& path = ridgeLines[i];
            if (path.size() > 0) {
                std::vector<QPointF> spath;
                spath.reserve(path.size());
                for (const Index2& p : path) {
                    spath.push_back(getScaledPoint(p, scaleLayerRidges));
                }
                painter.drawPolyline(&spath[0], int(spath.size()));
            }
        }

        // draw peaks?
        if (ui->ridges_showPeaks->isChecked()) {

            int numPeaks = int(prunedDivTree->peaks().size());
            std::vector<double> peakSizes(numPeaks, peakSize);
            if (ui->ridges_showProm->isChecked()) {
                std::vector<double> peakProm = terrainAnalysis->getPeaksProminence();
                for (int i = 0; i < numPeaks; i++) {
                    peakSizes[i] *= Math::Clamp(std::sqrt(peakProm[i]/100), 0.5, 4.0);
                }
            }
            painter.setPen(QPen(QBrush(QColor(0, 0, 0)), 2));
            painter.setBrush(QBrush(QColor(255, 128, 0)));
            for (unsigned int i = 0; i < prunedDivTree->peaks().size(); i++) {
                const Peak& p = prunedDivTree->peaks()[i];
                painter.drawEllipse(getScaledPoint(p.location, scaleLayerRidges), peakSizes[i], peakSizes[i]);
            }

            if (ui->ridges_displayPeakId->isChecked()) {
                painter.save();
                painter.setFont(QFont("Helvetica", 12, QFont::Bold));
                painter.setPen(QPen(QBrush(QColor(0, 0, 0)), 4));
                painter.scale(1, -1);
                for (unsigned int i = 0; i < prunedDivTree->peaks().size(); i++) {
                    const Peak& p = prunedDivTree->peaks()[i];
                    const QPointF pp = getScaledPoint(p.location, scaleLayerRidges);
                    painter.drawText(pp.x(), -pp.y(), "Peak " + QString::number(i)
                                     + ": " + QString::number(p.elevation));
                }
                painter.restore();
            }

            if (ui->ridges_showIsol->isChecked()) {
                std::vector<IsolationRecord> peakIsolation = terrainAnalysis->getPeaksIsolation();
                painter.setBrush(QBrush(QColor(64, 192, 64, 128)));
                painter.setPen(QPen(QBrush(QColor(0, 128, 0)), 4));

                int peakId = ui->ridges_showIsol_peak->value();
                if (peakId < int(peakIsolation.size())) {
                    if (peakIsolation[peakId].foundHigherGround) {
                        const Peak& p = prunedDivTree->peaks()[peakId];
                        QPointF pp = getScaledPoint(p.location, scaleLayerRidges);
                        QPointF pi = getScaledPoint(peakIsolation[peakId].closestHigherGround, scaleLayerRidges);
                        double isolSize = sqrt(double((pp.x() - pi.x()) * (pp.x() - pi.x()) + (pp.y() - pi.y()) * (pp.y() - pi.y())));
                        painter.drawEllipse(pp, isolSize, isolSize);
                        painter.drawLine(pp, pi);
                    }
                    else {
                        painter.drawRect(0, 0, scaleLayerRidges*hf.getSizeX(), scaleLayerRidges*hf.getSizeY());
                    }
                }
            }
        }
        // draw saddles?
        if (ui->ridges_showSaddles->isChecked())
        {
            painter.setPen(QPen(QBrush(QColor(0, 0, 0)), 2));
            painter.setBrush(QBrush(QColor(64, 192, 0)));
            for (const Saddle& s : prunedDivTree->saddles()) {
                if (!ignoreSea || s.elevation > seaLevel)
                    painter.drawEllipse(getScaledPoint(s.location, scaleLayerRidges), colSize, colSize);
            }
        }
    }
}
