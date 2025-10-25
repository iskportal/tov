plugins {
    kotlin("jvm") version "1.9.0"
    application
}





application {
    mainClass.set("tovextravaganza.TOVSolverKt")
}

repositories {
    mavenCentral()
}

dependencies {
    implementation(kotlin("stdlib"))
    implementation("commons-cli:commons-cli:1.6.0")
    implementation("org.apache.commons:commons-math3:3.6.1")
}

tasks.jar {
    manifest {
        attributes["Main-Class"] = "tovextravaganza.TOVSolverKt"
    }
    from(configurations.runtimeClasspath.get().map { if (it.isDirectory) it else zipTree(it) }) {
        exclude("META-INF/*.SF", "META-INF/*.DSA", "META-INF/*.RSA")
        duplicatesStrategy = DuplicatesStrategy.EXCLUDE
    }
}
