#include "terrainwidget.h"
#include <QMatrix4x4>
#include <QVector3D>
#include <QApplication>
#include <QCursor>
#include <iostream>


TerrainWidget::TerrainWidget(QWidget* parent) : QOpenGLWidget(parent)
{
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);

    cursorEnabled = false;
    cursorRadius = 1;
    cursorColor = Vector3(0,0,0);
}

TerrainWidget::~TerrainWidget()
{
    if (texId != 0) glDeleteTextures(1, &texId);
    if (skyboxVAO != 0) glDeleteVertexArrays(1, &skyboxVAO);
    if (bufferVerts != 0) glDeleteBuffers(1, &bufferVerts);
    if (bufferIndices != 0) glDeleteBuffers(1, &bufferIndices);
    if (meshVAO != 0) glDeleteVertexArrays(1, &meshVAO);
}

void TerrainWidget::initializeGL()
{
    initializeOpenGLFunctions();

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    shaderTerrain.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/shaders/terrain.vert");
    shaderTerrain.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/shaders/terrain.frag");
    shaderTerrain.link();
    shaderSkybox.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/shaders/skybox.vert");
    shaderSkybox.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/shaders/skybox.frag");
    shaderSkybox.link();

    glGenVertexArrays(1, &skyboxVAO);

    unsigned char foo = 255;
    glGenTextures(1, &texId);
    glBindTexture(GL_TEXTURE_2D, texId);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, (void*)(&foo));
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);

    emit glInitialized();
}

void TerrainWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, (GLint)w, (GLint)h);
}

void TerrainWidget::paintGL()
{
    double camDist = Norm(camera.getEye() - terrainBBox.center());
    double tradius = terrainBBox.radius();
    if (camDist < terrainBBox.radius())
        camera.setPlanes(100, 2*tradius);
    else
        camera.setPlanes(camDist - tradius, 2 * tradius + camDist);

    // Clear
    glClearColor(0.62f, 0.74f, 0.85f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    GLuint shader;

    // Sky
    shader = shaderSkybox.programId();
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glUseProgram(shader);
    glBindVertexArray(skyboxVAO);
    glUniform3f(glGetUniformLocation(shader, "CamPos"), camera.getEye()[0], camera.getEye()[1], camera.getEye()[2]);
    glUniform3f(glGetUniformLocation(shader, "CamLookAt"), camera.getAt()[0], camera.getAt()[1], camera.getAt()[2]);
    glUniform3f(glGetUniformLocation(shader, "CamUp"), camera.getUp()[0], camera.getUp()[1], camera.getUp()[2]);
    glUniform2f(glGetUniformLocation(shader, "iResolution"), width(), height());
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    // Terrain
	if (meshVAO > 0) {
        shader = shaderTerrain.programId();
        glUseProgram(shader);

        QMatrix4x4 matPerspective;
        matPerspective.perspective(Math::RadianToDegree(camera.getAngleOfViewV(width(), height())),
                        (GLdouble)width() / (GLdouble)height(),
                        camera.getNearPlane(), camera.getFarPlane());
        glUniformMatrix4fv(glGetUniformLocation(shader, "ProjectionMatrix"), 1, GL_FALSE, matPerspective.data());

        QMatrix4x4 matView;
        matView.lookAt(QVector3D(camera.getEye()[0], camera.getEye()[1], camera.getEye()[2]),
                       QVector3D(camera.getAt()[0],  camera.getAt()[1],  camera.getAt()[2]),
                       QVector3D(camera.getUp()[0],  camera.getUp()[1],  camera.getUp()[2]));
        glUniformMatrix4fv(glGetUniformLocation(shader, "ModelViewMatrix"), 1, GL_FALSE, matView.data());

		glUniform2f(glGetUniformLocation(shader, "u_worldMin"), terrainBBox.getMin()[0], terrainBBox.getMin()[1]);
        glUniform2f(glGetUniformLocation(shader, "u_worldSize"), terrainBBox.width(), terrainBBox.height());
        if (cursorEnabled) {
            glUniform1f(glGetUniformLocation(shader, "u_cursorRadius"), cursorRadius);
            glUniform4f(glGetUniformLocation(shader, "u_cursorColor"), cursorColor[0], cursorColor[1],
                                                                       cursorColor[2], 0.5);
            glUniform2f(glGetUniformLocation(shader, "u_cursorWorld"), cursorPos[0], cursorPos[1]);
        }
        else {
            glUniform1f(glGetUniformLocation(shader, "u_cursorRadius"), -1);
        }

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texId);
		glUniform1i(glGetUniformLocation(shader, "u_texture"), 0);

		glBindVertexArray(meshVAO);
		glDrawElements(GL_TRIANGLES, numTriangles*3, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);

		glUseProgram(0);
    }
}

void TerrainWidget::checkGLError()
{
    GLenum err;
    while ((err = glGetError()) != GL_NO_ERROR) {
        std::cerr << err << std::endl;
    }
}


void TerrainWidget::setHeightfield(const HeightField& hf, bool centerCamera)
{
	makeCurrent();

	int nx = hf.getSizeX();
	int ny = hf.getSizeY();

	// bbox
	Box2 bbox = hf.getDomain();
	double a, b;
	hf.getRange(a, b);
	terrainBBox = Box3(
		Vector3(bbox.getMin()[0], bbox.getMin()[1], a),
		Vector3(bbox.getMax()[0], bbox.getMax()[1], b));

	// camera
	if (centerCamera) camera = Camera::View(terrainBBox);

	// verts
	int idx = 0;
	std::vector<GLfloat> verts(3*nx*ny);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
            Vector3 p = hf.Vertex(i, j);
			verts[idx++] = GLfloat(p[0]);
			verts[idx++] = GLfloat(p[1]);
			verts[idx++] = GLfloat(p[2]);
		}
	}

	// tris
	numTriangles = (nx - 1) * (ny - 1) * 2;
	idx = 0;
	std::vector<GLuint> indices(numTriangles * 3);
	for (int i = 1; i < nx; i++) {
		for (int j = 1; j < ny; j++) {
			GLuint v00 = (i - 1) * ny + j - 1;
			GLuint v01 = (i - 1) * ny + j;
			GLuint v10 = i * ny + j - 1;
			GLuint v11 = i * ny + j;

			indices[idx++] = v00;
			indices[idx++] = v01;
			indices[idx++] = v10;

			indices[idx++] = v10;
			indices[idx++] = v01;
			indices[idx++] = v11;
		}
	}

	// update buffers
	if (bufferVerts > 0) glDeleteBuffers(1, &bufferVerts);
	if (bufferIndices > 0) glDeleteBuffers(1, &bufferIndices);
	if (meshVAO > 0) glDeleteVertexArrays(1, &meshVAO);

	glGenVertexArrays(1, &meshVAO);
	glBindVertexArray(meshVAO);

    glUseProgram(shaderTerrain.programId());
    GLuint attribVertexLoc = glGetAttribLocation(shaderTerrain.programId(), "a_position");

	glGenBuffers(1, &bufferVerts);
	glBindBuffer(GL_ARRAY_BUFFER, bufferVerts);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * verts.size(), &verts[0], GL_STATIC_DRAW);
	glVertexAttribPointer(attribVertexLoc, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(attribVertexLoc);

	glGenBuffers(1, &bufferIndices);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferIndices);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * indices.size(), &indices[0], GL_STATIC_DRAW);

	glBindVertexArray(0);
	glUseProgram(0);
}

void TerrainWidget::resetCamera()
{
    camera = Camera::View(terrainBBox);
    update();
}

void TerrainWidget::setTexture(const QImage& img)
{
	makeCurrent();

    glBindTexture(GL_TEXTURE_2D, texId);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width(), img.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, img.bits());
    glGenerateMipmap(GL_TEXTURE_2D);
}

void TerrainWidget::setCursorStyle(float radius, const Vector3 &color)
{
    cursorRadius = radius;
    cursorColor = color;
}

void TerrainWidget::showCursor(const Vector3 &worldPos)
{
    cursorEnabled = true;
    cursorPos = worldPos;
}

void TerrainWidget::hideCursor()
{
    cursorEnabled = false;
}


void TerrainWidget::mousePressEvent(QMouseEvent * e)
{
    x0 = e->globalPosition().x();
    y0 = e->globalPosition().y();

    if (e->buttons() & Qt::RightButton) {
        Ray ray = camera.pixelToRay(e->pos().x(), e->pos().y(), width(), height());
        emit rayQueried(ray);
    }

    update();
}

void TerrainWidget::mouseMoveEvent(QMouseEvent * e)
{
    int x = e->globalPosition().x();
    int y = e->globalPosition().y();

    double moveScale = Norm(camera.getEye() - camera.getAt()) * 0.015 * 0.05;

	if (e->buttons() & Qt::LeftButton) {
        if (e->modifiers() & Qt::ShiftModifier) {
            camera.backForth((y - y0) * moveScale);
            QApplication::setOverrideCursor(QCursor(Qt::SplitVCursor));
        }
        else if (e->modifiers() & Qt::AltModifier) {
            camera.leftRightPlane((x - x0) * moveScale);
            camera.upDownPlane((y - y0) * moveScale);
            QApplication::setOverrideCursor(QCursor(Qt::SizeAllCursor));
        }
        else {
            camera.leftRightRound((x0 - x) * 0.01);
            camera.upDownRound((y0 - y) * 0.005);
        }
	}
    else if (e->buttons() & Qt::MiddleButton) {
        camera.leftRightPlane((x - x0) * moveScale);
        camera.upDownPlane((y - y0) * moveScale);
        QApplication::setOverrideCursor(QCursor(Qt::SizeAllCursor));
    }
    else if (e->buttons() & Qt::RightButton) {
        Ray ray = camera.pixelToRay(e->pos().x(), e->pos().y(), width(), height());
        emit rayQueried(ray);
    }

    x0 = e->globalPosition().x();
    y0 = e->globalPosition().y();

    update();
}

void TerrainWidget::mouseReleaseEvent(QMouseEvent*)
{
    QApplication::setOverrideCursor(QCursor(Qt::ArrowCursor));    
	update();
}

void TerrainWidget::wheelEvent(QWheelEvent* e)
{
    int numDegrees = e->angleDelta().y() / 8;
    int numSteps = numDegrees / 15;
	
    if (!(e->modifiers()&Qt::ShiftModifier)) {
        Vector3 direction = camera.getAt() - camera.getEye();
        double currentDist = Norm(direction);
        direction = Normalized(direction);
        double speed = 0.1*currentDist;

        camera.setEye(camera.getEye() + direction * speed * numSteps);
    }

    update();
}
